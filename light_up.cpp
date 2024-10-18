#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>

using namespace std;

int len = 1000, init_len = 1000;
const int repeat_times_bound = 1000;
int repeat_times = 0;
double sum = 0, sum2 = 0;
const double T_0 = 100, cooling_rate = 1 - 1e-4, T_min = 0.01; 
double loss = 100000, loss_num = 100000;
int num = 0, sum_num = 0;
const double p_max = 6, p_min = 6, p_step = 1;
double p = 0, temp_p = 0, init_p = 0;
bool is_init = false;  //change this to initialize the p
bool is_stable = false;

struct vertex
{
    int state;
    vector<int> neighbor;
    vertex(): state(rand()%2), neighbor(vector<int>()) {}
};

struct edge
{
    bool is_exist = false;
    int value = false;
};

struct graph
{
    vector<vertex> vertices;
    vector<pair<int, int>> edges;
    vector<vector<edge>> edge_matrix;
    double loss_1, loss_21, loss_31;

    void add_edge(int u, int v)
    {
        vertices[u].neighbor.push_back(v);
        vertices[v].neighbor.push_back(u);
        edges.push_back(make_pair(u, v));
        edge_matrix[u][v].is_exist = true;
        edge_matrix[v][u].is_exist = true;
    }

    graph(int n): vertices(vector<vertex>(n))
    {
        edge_matrix = vector<vector<edge>>(n, vector<edge>(n));
        for (int i = 0; i < n; ++i)
        {
            for (int j = i + 1; j < n; ++j)
            {
                if (((double)rand() / RAND_MAX)<0.01)
                {
                    add_edge(i, j);
                }
            }
        }
        loss_1 = 0;
        loss_21 = 0;
        loss_31 = 0;
        for (vertex v: vertices)
        {
            loss_1 += v.state;
        }
        for (pair<int, int> e: edges)
        {
            edge_matrix[e.first][e.second].value = vertices[e.first].state^vertices[e.second].state;
            edge_matrix[e.second][e.first].value = vertices[e.first].state^vertices[e.second].state;
            loss_21 += double(edge_matrix[e.first][e.second].value);
            loss_31 += double(edge_matrix[e.first][e.second].value) * (1.0 / vertices[e.first].neighbor.size() + 1.0 / vertices[e.second].neighbor.size());
        }
    }

    void light_up(int u)
    {
        vertices[u].state = 1 - vertices[u].state;
        loss_1 += 2 * vertices[u].state - 1;
        for (int v: vertices[u].neighbor)
        {
            vertices[v].state = 1 - vertices[v].state;
            loss_1 += 2 * vertices[v].state - 1;
        }
        for (int v: vertices[u].neighbor)
            for (int w: vertices[v].neighbor)
            {
                if (w != u)
                {
                    if (edge_matrix[v][w].value != (vertices[v].state^vertices[w].state))
                    {
                        edge_matrix[v][w].value = 1 - edge_matrix[v][w].value;
                        edge_matrix[w][v].value = 1 - edge_matrix[w][v].value;
                        loss_21 += 2 * (edge_matrix[v][w].value) - 1;
                        loss_31 += (2 * (edge_matrix[v][w].value) - 1) * (1.0 / vertices[v].neighbor.size() + 1.0 / vertices[w].neighbor.size());
                    }
                }
            }
    }

    void print()
    {
        for (vertex v: vertices)
        {
            cout << v.state << " ";
        }
        cout << endl;
    }
}g(0);

void gradient_descent()
{
    bool is_converged = false;
    while (!is_converged)
    {
        is_converged = true;
        for (int i = 0; i < len; ++i)
        {
            sum_num++;
            g.light_up(i);
            double temp_loss_num = g.loss_1;
            if (temp_loss_num < loss_num)
            {
                loss_num = temp_loss_num;
                is_converged = false;
            }
            else
            {
                g.light_up(i);
            }
        }
    }
}

void main_algorithm()
{
    g = graph(len);
    loss = 100000;
    num = 0;
    sum_num = 0;
    temp_p = p;
    double T = T_0;
    while (true)
    {
        int x = rand() % len;
        g.light_up(x);
        double temp_loss = g.loss_1 + temp_p * g.loss_31; // change this to change the loss function
        if (loss >= temp_loss)
        {
            loss = temp_loss;
            num = 0;
        }
        else if (exp((loss - temp_loss) / T) > ((double)rand() / RAND_MAX))
        {
            loss = temp_loss;
            num = 0;
        }
        else
        {
            g.light_up(x);
        }
        num++;
        sum_num++;
        T *= cooling_rate;
        if (T < T_min)
        {
            break;
        }
        if (T < 1 && T > 0.1)
        {
            temp_p *= 1-2e-5;
        }
        if (T < 0.1)
        {
            temp_p *= 1-1e-4;
        }
    }
    loss_num = g.loss_1;
    gradient_descent();
    loss = g.loss_1 + temp_p * g.loss_31; // change this to change the loss function
    // print_light();
    cout << "len: " << len << "\tp: " << p << "\tloss_num: " << loss_num << " \tloss: " << loss << " \tsum_num: " << sum_num << endl;
    sum += loss_num;
    sum2 += loss_num * loss_num;
}

void repeat_algorithm()
{
    sum = 0;
    sum2 = 0;
    is_stable = false;
    for (repeat_times = 1; repeat_times <= repeat_times_bound; ++repeat_times)
    {
        main_algorithm();
        if (repeat_times >= 10)
            if (sum2 / (repeat_times*(repeat_times-1)) - ((sum*sum) / (repeat_times*repeat_times*(repeat_times-1))) < (sum / repeat_times * 0.005) * (sum / repeat_times * 0.005))
            {
                is_stable = true;
                break;
            }
    }
}

void extract_last_entry(int& len_value, double& p_value) {
    ifstream infile("light_output.txt");
    if (!infile) {
        cerr << "无法打开文件 light_output.txt" << endl;
        return;
    }

    string line;
    string last_line;
    while (getline(infile, line)) {
        if (!line.empty()) {
            last_line = line;
        }
    }

    infile.close();

    if (last_line.empty()) {
        cerr << "文件为空或没有有效数据" << endl;
        return;
    }

    istringstream iss(last_line);
    string part;
    while (getline(iss, part, '\t')) {
        if (part.find("len:") != string::npos) {
            len_value = stod(part.substr(part.find(":") + 1));
        } else if (part.find("p:") != string::npos) {
            p_value = stod(part.substr(part.find(":") + 1));
        }
    }
}

int main()
{
    srand(time(0));

    extract_last_entry(init_len, init_p);

    ofstream outfile("light_output.txt", std::ios::app);
    if (!outfile) {
        cerr << "无法打开文件 light_output.txt" << endl;
        return 1;
    }
    
    for (len = 1000; len<=1000; len+=10)
    {
        for (p = p_min; p <= p_max+1e-6; p += p_step)
        {
            if (is_init)
            {
                p = init_p+p_step;
                is_init = false;
                if (p>p_max)
                    break;
            }
            repeat_algorithm();
            outfile << "len: " << len << "\tp: " << p << "\taverage loss_num: " << sum / repeat_times;
            if (is_stable)
            {
                outfile << endl;
            }
            else
            {
                outfile << "\tnot stable: " << (sqrt(sum2 / (repeat_times*(repeat_times-1)) - ((sum*sum) / (repeat_times*repeat_times*(repeat_times-1))))/(sum / repeat_times)) << endl;
            }
        }
    }

    // outfile << "这是写入文件的第一行。" << endl;
    // outfile << "这是写入文件的第二行。" << endl;

    outfile.close();

    return 0;
}
