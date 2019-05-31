#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<fstream>
#include<iostream>
#include<cstring>
#include<random>
#include<string>
#include<math.h>
#include<algorithm>
#include<queue>

using namespace std;

int dim;//总的节点数目

int node_num;//必须经过节点的数目

int edge_num;//必须经过边的条数

const int infinite = 100;//邻接表中表示两节点不直接相连

int** Adjacency_Matrix;//邻接矩阵

int* node;//必须经过节点数组,第一个元素是起点，最后一个元素是终点

vector<int> edge;//必须经过的边数组，外层必须经过边的条数，内层下标0起点，下标1终点

int** W_Matrix; //Floyd算法所需邻接表

int** R_Matrix;//Floyd算法所需路由表

int** A_Matrix;//在最短路径基础上只有必经节点的邻接表

vector<int>** Path;//在最短路径基础上必经节点之间的最短路径

int*** Dp;//状态矩阵，外层下标代表已经经过节点状态，中层代表当前状态下到达的最后一个节点，内层代表是否经过必经边，数值为当前的最短路径

vector<int>*** State_Path;//当前状态的最短路径

int ans = -1;//最后得出的最短路径值

vector<int> ans_path;//最后得出的最短路径

int max_hop;//允许经过的最大节点数目

int min_hop = -1;//计算得到的满足除节点数目之外所有条件的最小跳数

int shortest_path_hop = -1;//计算得到的满足除节点数目之外的所有条件的最短路径

vector<vector<int>>** possible_path;//外两层下标为节点号，内层外层vector为可能的路径，内层内层vector为路径经过节点

vector<int>** possible_path_length;//对应于上面每个路径的长度

int**** CDp;//加上计算节点数的状态矩阵，前三层与Dp相同，最内层为跳数

vector<int>**** CState_Path;//同State_Path

struct Node{
	int id;
	int layer;
	int path_length;
	vector<int> path;
	Node(int a, int b, int c):id(a),layer(b),path_length(c){}
	Node() {
		id = -1; layer = -1; path_length = 100;
	}
};

void Iniinitialize() {
	//打开输入文件
	cout << "输入文件路径:";
	string file_path;
	cin >> file_path;

	ifstream input(file_path, ios::in);

	if (!input) { cout << "输入文件出错"; }

	//建立邻接表
	input >> dim;//输入总结点数目

	Adjacency_Matrix = new int*[dim];
	for (int i = 0; i < dim; i++) {
		Adjacency_Matrix[i] = new int[dim];
	}

	//邻接表赋值
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			input >> Adjacency_Matrix[i][j];
		}
	}

	for (int i = 0; i < dim; i++) {
		for (int j = i; j < dim; j++) {
			if (i == j) {
				if (Adjacency_Matrix[i][j] != 0) throw logic_error("error");
			}
			else {
				if (Adjacency_Matrix[i][j] != Adjacency_Matrix[j][i]) throw logic_error("error");
			}
		}
	}

	//建立必须经过节点的数组
	input >> node_num;//输入必须经过节点的数目
	node = new int[node_num];

	//起点第一个，终点最后一个，其余无所谓
	for (int i = 0; i < node_num; i++) {
		input >> node[i];
	}

	//建立必须经过的边的数组
	input >> edge_num;
	int tempt_node_id;
	for (int i = 0; i < edge_num * 2; i++) {
		input >> tempt_node_id;
		edge.push_back(tempt_node_id);
	}

	//输入要求经过最多节点数
	input >> max_hop;
}

void Floyd(){
	//复制邻接表
	W_Matrix = new int*[dim];
	for (int i = 0; i < dim; i++) {
		W_Matrix[i] = new int[dim];
	}

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			W_Matrix[i][j] = Adjacency_Matrix[i][j];
		}
	}

	//建立路由表
	R_Matrix = new int*[dim];
	for (int i = 0; i < dim; i++) {
		R_Matrix[i] = new int[dim];
	}

	//初始化路由表
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
		{
			R_Matrix[i][j] = j;
		}

	//算法迭代
	for (int k = 0; k < dim; k++)
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				if (W_Matrix[i][j] > W_Matrix[i][k] + W_Matrix[k][j])
				{
					W_Matrix[i][j] = W_Matrix[i][k] + W_Matrix[k][j];
					R_Matrix[i][j] = R_Matrix[i][k];
				}
}

void Format() {
	//建立只有必经点的新的邻接矩阵
	A_Matrix = new int*[node_num];
	for (int i = 0; i < node_num; i++) {
		A_Matrix[i] = new int[node_num];
	}

	for (int i = 0; i < node_num; i++) {
		for (int j = 0; j < node_num; j++) {
			A_Matrix[i][j] = W_Matrix[node[i]][node[j]];
		}
	}


	//记录新的邻接矩阵中两点间的路径
	Path = new vector<int>*[node_num];
	for (int i = 0; i < node_num; i++) {
		Path[i] = new vector<int>[node_num];
		for (int j = 0; j < node_num; j++) {
			int temp = R_Matrix[node[i]][node[j]];
			while (temp != node[j]) {
				Path[i][j].push_back(temp);
				temp = R_Matrix[temp][node[j]];
			}
			Path[i][j].push_back(node[j]);
		}
	}

	//将必经边置入
	for (int i = 0; i < edge_num*2; i=i+2) {
		int start_node = -1;
		int end_node = -1;
		for (int j = 0; j < node_num; j++) {
			if (node[j] == edge[i]) start_node = j;
			if (node[j] == edge[i + 1]) end_node = j;
		}

		A_Matrix[start_node][end_node] = Adjacency_Matrix[node[start_node]][node[end_node]];
		A_Matrix[end_node][start_node] = Adjacency_Matrix[node[end_node]][node[start_node]];
		Path[start_node][end_node].clear();
		Path[end_node][start_node].clear();
		Path[start_node][end_node].push_back(node[end_node]);
		Path[end_node][start_node].push_back(node[start_node]);

	}
}

void Print() {
	{
		cout << "Adjacency_Matrix:" << endl;
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				cout << Adjacency_Matrix[i][j] << " ";
			}
			cout << endl;
		}

		cout << "node:" << endl;
		for (int i = 0; i < node_num; i++) {
			cout << node[i] << " ";
		}
		cout << endl;

		cout << "W_Matrix:" << endl;
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				cout << W_Matrix[i][j] << " ";
			}
			cout << endl;
		}

		cout << "R_Matrix:" << endl;
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				cout << R_Matrix[i][j] << " ";
			}
			cout << endl;
		}

		cout << "A_Matrix:" << endl;
		for (int i = 0; i < node_num; i++) {
			for (int j = 0; j < node_num; j++) {
				cout << A_Matrix[i][j] << " ";
			}
			cout << endl;
		}

		cout << "edge:" << endl;
		for (int i = 0; i < edge_num * 2; i = i + 2) {
			cout << edge[i] << "-" << edge[i + 1] << " ";
		}
	}
}

int Min(int a, int b) {
	if (a == -1) return b;
	if (b == -1) return a;
	return a < b ? a : b;

}

void DP() {
	int state_num = pow(2, node_num);//状态数

	int sub_state_num = pow(2, edge_num);//子状态数，主要取决于必经点的个数

	//建立状态矩阵
	Dp = new int**[state_num];
	for (int i = 0; i < state_num; i++) {
		Dp[i] = new int*[node_num];
		for (int j = 0; j < node_num; j++) {
			Dp[i][j] = new int[sub_state_num];
			for (int k = 0; k < sub_state_num; k++) {
				Dp[i][j][k] = -1;
			}
		}
	}

	Dp[1 << 0][0][0] = 0;//选择起点

	State_Path = new vector<int>**[state_num];
	for (int i = 0; i < state_num; i++) {
		State_Path[i] = new vector<int>*[node_num];
		for (int j = 0; j < node_num; j++) {
			State_Path[i][j] = new vector<int>[sub_state_num];
		}
	}

	State_Path[1 << 0][0][0].push_back(0);//初始状态路径加入起点

	//DP状态压缩算法
	for (int i = 0; i < 1 << node_num; i++)
	{
		for (int j = 0; j < node_num; j++)
		{
			for (int l = 0; l < 1 << edge_num; l++) {
				if (Dp[i][j][l] != -1)
				{
					for (int k = 0; k < node_num; k++)
					{
						//相同状态更新最短路径
						if (Dp[i | (1 << k)][k][l] == -1 || Dp[i | (1 << k)][k][l] > (Dp[i][j][l] + A_Matrix[j][k])) {
							Dp[i | (1 << k)][k][l] = Dp[i][j][l] + A_Matrix[j][k];
							State_Path[i | (1 << k)][k][l].assign(State_Path[i][j][l].begin(), State_Path[i][j][l].end());
							vector<int>::iterator it2 = Path[j][k].begin();
							while (it2 != Path[j][k].end()) {
								State_Path[i | (1 << k)][k][l].push_back(*it2);
								it2++;
							}
						}

						//如果状态发生改变，则除更新自己状态以外，尝试更新新的状态
						vector<int>::iterator it = edge.begin();
						vector<int>::iterator it1 = find(edge.begin(), edge.end(), node[j]);
						while (it1 != edge.end())
						{
							int another_node;
							if ((it1 - edge.begin()) % 2 == 0) another_node = *(it1 + 1);
							else another_node = *(it1 - 1);

							if (node[k] == another_node) {
								int new_sub_state = 1 << ((it1 - edge.begin()) / 2);
								if (Dp[i | (1 << k)][k][l | new_sub_state] == -1 || Dp[i | (1 << k)][k][l | new_sub_state] > (Dp[i][j][l] + A_Matrix[j][k])) {
									Dp[i | (1 << k)][k][l | new_sub_state] = Dp[i][j][l] + A_Matrix[j][k];
									State_Path[i | (1 << k)][k][l | new_sub_state].assign(State_Path[i][j][l].begin(), State_Path[i][j][l].end());
									vector<int>::iterator it3 = Path[j][k].begin();
									while (it3 != Path[j][k].end()) {
										State_Path[i | (1 << k)][k][l | new_sub_state].push_back(*it3);
										it3++;
									}
								}
							}

							it1 = find(it1 + 1, edge.end(), node[j]);
						}

						//将暂时路径最小切满足条件的状态计入结果
						if ((i | (1 << k)) == (1 << node_num) - 1 && k == node_num - 1 && l == sub_state_num-1) {
							if (ans == -1 || ans > Dp[i | (1 << k)][k][l]) {
								ans = Dp[i | (1 << k)][k][l];
								ans_path.assign(State_Path[i | (1 << k)][k][l].begin(), State_Path[i | (1 << k)][k][l].end());
							}
						}
					}
				}
			}
		}
	}
	shortest_path_hop = ans_path.size();
}

void Calculate_Min_Hop(){
	int** cW_Matrix= new int*[dim];
	for (int i = 0; i < dim; i++) {
		cW_Matrix[i] = new int[dim];
	}

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (Adjacency_Matrix[i][j] != 100 && i != j)
				cW_Matrix[i][j] = 1;
			else 
				cW_Matrix[i][j] = Adjacency_Matrix[i][j];
		}
	}

	for (int k = 0; k < dim; k++)
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				if (cW_Matrix[i][j] > cW_Matrix[i][k] + cW_Matrix[k][j])
				{
					cW_Matrix[i][j] = cW_Matrix[i][k] + cW_Matrix[k][j];
				}

	int** cA_Matrix = new int*[node_num];
	for (int i = 0; i < node_num; i++) {
		cA_Matrix[i] = new int[node_num];
	}

	for (int i = 0; i < node_num; i++) {
		for (int j = 0; j < node_num; j++) {
			cA_Matrix[i][j] = cW_Matrix[node[i]][node[j]];
		}
	}

	int state_num = pow(2, node_num);

	int sub_state_num = pow(2, edge_num);
					 
	int*** cDp = new int**[state_num];
	for (int i = 0; i < state_num; i++) {
		cDp[i] = new int*[node_num];
		for (int j = 0; j < node_num; j++) {
			cDp[i][j] = new int[sub_state_num];
			for (int k = 0; k < sub_state_num; k++) {
				cDp[i][j][k] = -1;
			}
		}
	}

	cDp[1 << 0][0][0] = 0;

	for (int i = 0; i < 1 << node_num; i++)
	{
		for (int j = 0; j < node_num; j++)
		{
			for (int l = 0; l < 1 << edge_num; l++) {
				if (cDp[i][j][l] != -1)
				{
					for (int k = 0; k < node_num; k++)
					{
						//相同状态更新最短路径
						if (cDp[i | (1 << k)][k][l] == -1 || cDp[i | (1 << k)][k][l] >(cDp[i][j][l] + cA_Matrix[j][k])) {
							cDp[i | (1 << k)][k][l] = cDp[i][j][l] + cA_Matrix[j][k];
						}

						//如果状态发生改变，则除更新自己状态以外，尝试更新新的状态
						vector<int>::iterator it = edge.begin();
						vector<int>::iterator it1 = find(edge.begin(), edge.end(), node[j]);
						while (it1 != edge.end())
						{
							int another_node;
							if ((it1 - edge.begin()) % 2 == 0) another_node = *(it1 + 1);
							else another_node = *(it1 - 1);

							if (node[k] == another_node) {
								int new_sub_state = 1 << ((it1 - edge.begin()) / 2);
								if (cDp[i | (1 << k)][k][l | new_sub_state] == -1 || cDp[i | (1 << k)][k][l | new_sub_state] > (cDp[i][j][l] + cA_Matrix[j][k])) {
									cDp[i | (1 << k)][k][l | new_sub_state] = cDp[i][j][l] + cA_Matrix[j][k];
								}
							}

							it1 = find(it1 + 1, edge.end(), node[j]);
						}

						//将暂时路径最小切满足条件的状态计入结果
						if ((i | (1 << k)) == (1 << node_num) - 1 && k == node_num - 1 && l == sub_state_num - 1) {
							min_hop = Min(min_hop, cDp[i | (1 << k)][k][l]);
						}
					}
				}
			}
		}
	}
	min_hop = min_hop + 1;
}

void BFS() {
	possible_path = new vector<vector<int>>*[node_num];
	for (int i = 0; i < node_num; i++) {
		possible_path[i] = new vector<vector<int>>[node_num];
	}

	possible_path_length = new vector<int>*[node_num];
	for (int i = 0; i < node_num; i++) {
		possible_path_length[i] = new vector<int>[node_num];
	}

	for (int i = 0; i < node_num; i++) {
		for (int j = 0; j < node_num; j++) {
			int s_path_hop = Path[i][j].size();
			if (s_path_hop == 1) {
				possible_path[i][j].push_back(Path[i][j]);
				possible_path_length[i][j].push_back(A_Matrix[i][j]);
				continue;
			}
				
			s_path_hop = s_path_hop - 1;
			Node* possible_node = new Node[s_path_hop];

			queue<Node> q;
			Node nod(node[i], 0, 0);
			q.push(nod);

			while (!q.empty()) {
				if (nod.layer > s_path_hop) throw logic_error("error");
				if (nod.layer > s_path_hop) break;
				nod = q.front();
				q.pop();

				for (int k = 0; k < dim; k++) {
					if (Adjacency_Matrix[nod.id][k] != 0 && Adjacency_Matrix[nod.id][k] != 100) {
						if (k == node[j]) {
							if (nod.path_length < possible_node[nod.layer].path_length) {
								possible_node[nod.layer].id = nod.id;
								possible_node[nod.layer].layer = nod.layer;
								possible_node[nod.layer].path_length = nod.path_length + Adjacency_Matrix[nod.id][k];
								possible_node[nod.layer].path = nod.path;
								possible_node[nod.layer].path.push_back(node[j]);
							}
						}
						else if (nod.layer < s_path_hop - 1) {
							Node next(k, nod.layer + 1, nod.path_length + Adjacency_Matrix[nod.id][k]);
							next.path = nod.path;
							next.path.push_back(k);
							q.push(next);
						}
					}
				}
			}

			int current_min = 100;
			for (int l = 0; l < s_path_hop; l++) {
				if (possible_node[l].id != -1) {
					if (possible_node[l].path_length < current_min) {
						possible_path[i][j].push_back(possible_node[l].path);
						possible_path_length[i][j].push_back(possible_node[l].path_length);
						current_min = possible_node[l].path_length;
					}
				}
			}
			possible_path[i][j].push_back(Path[i][j]);
			possible_path_length[i][j].push_back(A_Matrix[i][j]);
		}
	}
}

void CDP() {
	int state_num = pow(2, node_num);//总状态

	int sub_state_num = pow(2, edge_num);//子状态，与必经边数量有关

	CDp = new int***[state_num];
	for (int i = 0; i < state_num; i++) {
		CDp[i] = new int**[node_num];
		for (int j = 0; j < node_num; j++) {
			CDp[i][j] = new int*[sub_state_num];
			for (int k = 0; k < sub_state_num; k++) {
				CDp[i][j][k] = new int[max_hop];
				for (int l = 0; l < max_hop; l++) {
					CDp[i][j][k][l] = -1;
				}
			}
		}
	}

	CDp[1 << 0][0][0][0] = 0;//选择起点

	CState_Path = new vector<int>***[state_num];
	for (int i = 0; i < state_num; i++) {
		CState_Path[i] = new vector<int>**[node_num];
		for (int j = 0; j < node_num; j++) {
			CState_Path[i][j] = new vector<int>*[sub_state_num];
			for (int k = 0; k < sub_state_num; k++) {
				CState_Path[i][j][k] = new vector<int>[max_hop];
			}
		}
	}

	CState_Path[1 << 0][0][0][0].push_back(0);//初始状态路径加入起点

	//DP状态压缩算法
	for (int i = 0; i < 1 << node_num; i++)
	{
		for (int j = 0; j < node_num; j++)
		{
			for (int l = 0; l < 1 << edge_num; l++) {
				for (int m = 0; m < max_hop; m++) {
					if (CDp[i][j][l][m] != -1)
					{
						for (int k = 0; k < node_num; k++)
						{
							for (int n = 0; n < possible_path[j][k].size(); n++) {
								if (m + possible_path[j][k][n].size() < max_hop) {
									//相同状态更新最短路径
									if (CDp[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ] == -1 || CDp[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ] > (CDp[i][j][l][m] + possible_path_length[j][k][n])) {
										CDp[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ] = (CDp[i][j][l][m] + possible_path_length[j][k][n]);
										CState_Path[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ].assign(CState_Path[i][j][l][m].begin(), CState_Path[i][j][l][m].end());
										vector<int>::iterator it2 = possible_path[j][k][n].begin();
										while (it2 != possible_path[j][k][n].end()) {
											CState_Path[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ].push_back(*it2);
											it2++;
										}
									}

									//如果状态发生改变，则除更新自己状态以外，尝试更新新的状态
									vector<int>::iterator it = edge.begin();
									vector<int>::iterator it1 = find(edge.begin(), edge.end(), node[j]);
									while (it1 != edge.end())
									{
										int another_node;
										if ((it1 - edge.begin()) % 2 == 0) another_node = *(it1 + 1);
										else another_node = *(it1 - 1);

										if (node[k] == another_node) {
											int new_sub_state = 1 << ((it1 - edge.begin()) / 2);
											if (CDp[i | (1 << k)][k][l | new_sub_state][m + possible_path[j][k][n].size() ] == -1 || CDp[i | (1 << k)][k][l | new_sub_state][m + possible_path[j][k][n].size() ] > (CDp[i][j][l][m] + possible_path_length[j][k][n])) {
												CDp[i | (1 << k)][k][l | new_sub_state][m + possible_path[j][k][n].size() ] = CDp[i][j][l][m] + possible_path_length[j][k][n];
												CState_Path[i | (1 << k)][k][l | new_sub_state][m + possible_path[j][k][n].size() ].assign(CState_Path[i][j][l][m].begin(), CState_Path[i][j][l][m].end());
												vector<int>::iterator it3 = possible_path[j][k][n].begin();
												while (it3 != possible_path[j][k][n].end()) {
													CState_Path[i | (1 << k)][k][l | new_sub_state][m + possible_path[j][k][n].size() ].push_back(*it3);
													it3++;
												}
											}
										}

										it1 = find(it1 + 1, edge.end(), node[j]);
									}

									//将暂时路径最小切满足条件的状态计入结果
									if ((i | (1 << k)) == (1 << node_num) - 1 && k == node_num - 1 && l == sub_state_num - 1 && m + possible_path[j][k][n].size() < max_hop) {
										if (ans == -1 || ans > CDp[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ]) {
											ans = CDp[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ];
											ans_path.assign(CState_Path[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ].begin(), CState_Path[i | (1 << k)][k][l][m + possible_path[j][k][n].size() ].end());
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void main() {
	Iniinitialize();//初始化

	Calculate_Min_Hop();//计算满足除了节点数目以外的所有条件的经过最少节点数

	Floyd();//弗洛伊德算法计算必经节点间的最短距离

	Format();//以必经节点为元素重新构造邻接矩阵

	DP();//状态压缩算法计算最短路径

	cout << "总结点数：" << dim << endl;
	cout << "必经节点：";
	for (int i = 0; i < node_num; i++) {
		cout << node[i] << " ";
	}
	cout << endl;
	cout << "必经边：";
	for (int i = 0; i < edge_num * 2; i=i+2) {
		cout << edge[i] << "->" << edge[i + 1] << " ";
	}
	cout << endl;
	cout << "经过最大节点数目：" << max_hop;
	cout << endl << "************************************************************" << endl;

	if (min_hop > max_hop){
		cout << "无解，输出次优解（即只不满足跳数限制其余条件均满足）";
		cout << endl << "最短路径：";
		vector<int>::iterator ite = ans_path.begin();
		while (ite != ans_path.end()) {
			if(ite!=ans_path.end()-1)
			cout << *ite << "->";
			else cout << *ite << endl;
			ite++;
		}
		cout << "最短路径长度：" << ans << endl;
	}

	else if (shortest_path_hop <= max_hop) {
		cout << "有解:";
		cout << endl << "最短路径：";
		vector<int>::iterator ite = ans_path.begin();
		while (ite != ans_path.end()) {
			if (ite != ans_path.end() - 1)
				cout << *ite << "->";
			else cout << *ite << endl;
			ite++;
		}
		cout << "最短路径长度：" << ans << endl;
	}

	else {
		ans = -1;
		ans_path.clear();
		BFS();
		CDP();
		cout << "有解:";
		cout << endl << "最短路径：";
		vector<int>::iterator iter = ans_path.begin();
		while (iter != ans_path.end()) {
			if (iter != ans_path.end() - 1)
				cout << *iter << "->";
			else cout << *iter << endl;
			iter++;
		}
		cout << "最短路径长度：" << ans << endl;
	}
	cout << endl;
	system("pause");
}