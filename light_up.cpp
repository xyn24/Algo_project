#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <utility>
#include <chrono>
#include <iomanip>
#include <random>

#define MAXN 1024
#define MAXEDGES 32768
#define MAXNB 128

using namespace std;

int len = 1000;
int repeat_times = 0;
const int repeat_times_bound = 1000;
double sum = 0, sum2 = 0, sum_time = 0;
const double T_0 = 5, cooling_rate = 1e-5, T_min = 0.01;
double T = 0;
const double relative_uncertainty = 0.005, edge_generate_probability = (1.0/5000);
long double loss = 100000, loss_num = 100000, temp_loss = 0, min_loss = 100000;
int num = 0, sum_num = 0;
int loss_type = 1;    // change this to change the loss function
bool is_stable = false;

random_device rd;                               // 用于生成种子
mt19937 gen(time(0));                           // 使用 Mersenne Twister 生成器
uniform_int_distribution<> dis(1, 100000);      // 定义一个均匀分布，范围为 1 到 100000
uniform_real_distribution<> dis_real(0.0, 1.0); // 定义一个均匀分布，范围为 0.0 到 1.0
uniform_int_distribution<> dis_bool(0, 1);      // 定义一个均匀分布，范围为 0 到 1

struct vertex
{
    int state;
    int neighbor[MAXNB];
    int size;
    int size_3;
    vertex() : state(dis_bool(gen)), size(0) {}

    vertex &operator=(const vertex &a)
    {
        state = a.state;
        size = a.size;
        size_3 = a.size_3;
        for (int i = 0; i < size; ++i)
        {
            neighbor[i] = a.neighbor[i];
        }
        return *this;
    }
};

struct graph
{
    vertex vertices[MAXN];
    int edges[MAXEDGES][2];
    int vertices_size, edge_size;
    bool edge_matrix[MAXN][MAXN];
    long double loss_1;

    void add_edge(int u, int v)
    {
        vertices[u].neighbor[vertices[u].size++] = v;
        vertices[v].neighbor[vertices[v].size++] = u;
        edges[edge_size][0] = u;
        edges[edge_size][1] = v;
        edge_size++;
        edge_matrix[u][v] = true;
        edge_matrix[v][u] = true;
    }

    graph(int n) : vertices_size(n), edge_size(0), loss_1(0)
    {
        for (int i = 0; i < n; ++i)
            vertices[i] = vertex();
        for (int i = 0; i < n; ++i)
            fill(edge_matrix[i], edge_matrix[i] + n, false);
        for (int i = 0; i < n; ++i)
        {
            for (int j = i + 1; j < n; ++j)
            {
                if ((dis_real(gen)) < edge_generate_probability)
                {
                    add_edge(i, j);
                }
            }
        }
        for (int i = 0; i < 20; ++i)
            for (int j = 0; j < 49; ++j)
                if (edge_matrix[50*i+j][50*i+j+1] == false)
                    add_edge(50*i+j, 50*i+j+1);
        for (int i = 0; i < 19; ++i)
            for (int j = 0; j < 50; ++j)
                if (edge_matrix[50*i+j][50*(i+1)+j] == false)
                    add_edge(50*i+j, 50*(i+1)+j);
        for (int i = 0; i < vertices_size; ++i)
        {
            vertices[i].size_3 = pow(vertices[i].size, 3);
            loss_1 += vertices[i].state;
        }
    }

    graph &operator=(const graph &a)
    {
        if (this == &a)
            return *this;
        vertices_size = a.vertices_size;
        edge_size = a.edge_size;
        loss_1 = a.loss_1;
        for (int i = 0; i < vertices_size; ++i)
        {
            vertices[i] = a.vertices[i];
        }
        for (int i = 0; i < edge_size; ++i)
        {
            edges[i][0] = a.edges[i][0];
            edges[i][1] = a.edges[i][1];
        }
        for (int i = 0; i < vertices_size; ++i)
        {
            for (int j = 0; j < vertices_size; ++j)
            {
                edge_matrix[i][j] = a.edge_matrix[i][j];
            }
        }
        return *this;
    }

    void light_up(int u)
    {
        vertices[u].state = 1 - vertices[u].state;
        loss_1 += 2 * vertices[u].state - 1;
        for (int i = 0; i < vertices[u].size; ++i)
        {
            int v = vertices[u].neighbor[i];
            vertices[v].state = 1 - vertices[v].state;
            loss_1 += 2 * vertices[v].state - 1;
        }
        for (int i = 0; i < vertices[u].size; ++i)
        {
            int v = vertices[u].neighbor[i];
            for (int j = 0; j < vertices[v].size; ++j)
            {
                int w = vertices[v].neighbor[j];
            }
        }
    }

    long double loss(int choose)
    {
        if (choose == 1)
            return loss_1;
        exit(2);
    }

    void print()
    {   
        for (int i = 0; i < vertices_size; ++i)
        {
            cout << vertices[i].state << " ";
        }
        cout << endl;
    }
} g(0), g_min(0);

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

void update()
{
    loss = temp_loss;
    // if (g.loss_1 < min_loss)
    // {
    //     min_loss = g.loss_1;
    //     g_min = g;
    // }
}

void main_algorithm()
{
    g = graph(len);
    auto start_time = chrono::high_resolution_clock::now();
    loss = 100000;
    min_loss = 100000;
    // loss = g.loss(loss_type);
    num = 0;
    sum_num = 0;
    for (int layer = 0; layer < 20; ++layer)
    {
        T = T_0;
        while (true)
        {
            loss = g.loss(loss_type);
            int x = dis(gen) % 50;
            x += 50 * layer;
            int d_loss = 0;
            for (int i = 0; i <= g.vertices[x].size; ++i)
            {
                int y;
                if (i == g.vertices[x].size)
                    y = x;
                else
                    y = g.vertices[x].neighbor[i];
                bool is_end = true;
                for (int j = 0; j < g.vertices[y].size; ++j)
                {
                    int z = g.vertices[y].neighbor[j];
                    if (z >= 50 * (layer + 1))
                    {
                        is_end = false;
                        break;
                    }
                }
                if (is_end)
                {
                    d_loss += 1 - 2 * g.vertices[y].state;
                }
            }
            if (0 >= d_loss)
            {
                // T *= 1 + 0.2 * (temp_loss - loss) / loss;
                update();
                g.light_up(x);
            }
            else if (exp(-d_loss / T) > (dis_real(gen)))
            {
                // T *= 1 + 0.1 * (temp_loss - loss) / loss;
                update();
                g.light_up(x);
            }
            num++;
            sum_num++;
            // T -= cooling_rate;
            // T *= cooling_rate;
            T *= 1 / (1 + cooling_rate * log(1 + sum_num));
            if (T < T_min)
            {
                break;
            }
            // if (sum_num % 20000 == 0)
            //     cout << "\tloss_num: " << g.loss_1 << "\tloss: " << loss << endl;
        }
    }
    // g = g_min;
    loss_num = g.loss_1;
    // gradient_descent();
    loss = g.loss(loss_type);
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> runtime = end_time - start_time;
    // g.print();
    cout << "len: " << len << "\tloss_num: " << loss_num << "\tloss: " << loss << " \tsum_num: " << sum_num << " \truntime: " << runtime.count() << " \tmin_loss: " << min_loss << endl;
    sum += loss_num;
    sum2 += loss_num * loss_num;
    sum_time += runtime.count();
}

void repeat_algorithm()
{
    sum = 0;
    sum2 = 0;
    sum_time = 0;
    is_stable = false;
    for (repeat_times = 1; repeat_times <= repeat_times_bound; ++repeat_times)
    {
        main_algorithm();
        if (repeat_times >= 10)
            if (sum2 / (repeat_times * (repeat_times - 1)) - ((sum * sum) / (repeat_times * repeat_times * (repeat_times - 1))) < (sum / repeat_times * relative_uncertainty) * (sum / repeat_times * relative_uncertainty))
            {
                is_stable = true;
                break;
            }
    }
}

void extract_last_entry(int &len_value, double &p_value)
{
    ifstream infile("light_output.txt");
    if (!infile)
    {
        cerr << "无法打开文件 light_output.txt" << endl;
        exit(1);
    }

    string line;
    string last_line;
    while (getline(infile, line))
    {
        if (!line.empty())
        {
            last_line = line;
        }
    }

    infile.close();

    if (last_line.empty())
    {
        cerr << "文件为空或没有有效数据" << endl;
        return;
    }

    istringstream iss(last_line);
    string part;
    while (getline(iss, part, '\t'))
    {
        if (part.find("len:") != string::npos)
        {
            len_value = stod(part.substr(part.find(":") + 1));
        }
        else if (part.find("p:") != string::npos)
        {
            p_value = stod(part.substr(part.find(":") + 1));
        }
    }
}

void main_function()
{
    // extract_last_entry(init_len, init_p);

    ofstream outfile("light_output.txt", std::ios::app);
    if (!outfile)
    {
        cerr << "无法打开文件 light_output.txt" << endl;
        exit(1);
    }

    repeat_algorithm();
    outfile << "len: " << len << "\taverage loss_num: " << sum / repeat_times << "\taverage runtime: " << sum_time / repeat_times;
    if (is_stable)
    {
        outfile << endl;
    }
    else
    {
        outfile << "\tnot stable: " << (sqrt(sum2 / (repeat_times * (repeat_times - 1)) - ((sum * sum) / (repeat_times * repeat_times * (repeat_times - 1)))) / (sum / repeat_times)) << endl;
    }

    outfile.close();
}

int main()
{
    main_function();

    return 0;
}