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

int len = 1000, init_len = 1000;
const int len_min = 1000, len_max = 1000, len_step = 10; //change this to change the len
int repeat_times = 0;
const int repeat_times_bound = 1000;
double sum = 0, sum2 = 0, sum_time = 0;
const double T_0 = 100, cooling_rate = 1 - 1e-4, T_min = 0.01;
double T = 0;
const double relative_uncertainty = 0.01, edge_generate_probability = 0.02;
long double loss = 100000, loss_num = 100000, temp_loss = 0;
int num = 0, sum_num = 0;
const double p_max = 1000, p_min = 60, p_step = 0.5; //change this to change the p
double p = 0, temp_p = 0, init_p = 0;
double p_momentum = 0, p_upbound = 100, p_lowbound = 0;
int loss_type = 4; // change this to change the loss function
bool is_init = false; // change this to initialize the p
bool is_stable = false;

random_device rd; // 用于生成种子
mt19937 gen(time(0)); // 使用 Mersenne Twister 生成器
uniform_int_distribution<> dis(1, 100000); // 定义一个均匀分布，范围为 1 到 100
uniform_real_distribution<> dis_real(0.0, 1.0); // 定义一个均匀分布，范围为 0.0 到 1.0
uniform_int_distribution<> dis_bool(0, 1); // 定义一个均匀分布，范围为 0 到 1

struct vertex
{
    int state;
    int neighbor[MAXNB];
    int size;
    vertex() : state(dis_bool(gen)), size(0) {}
};

struct graph
{
    vertex vertices[MAXN];
    int edges[MAXEDGES][2];
    int vertices_size, edge_size;
    bool edge_matrix[MAXN][MAXN];
    long double loss_1, loss_21, loss_31, loss_41;

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

    graph(int n) : vertices_size(n), edge_size(0), loss_1(0), loss_21(0), loss_31(0), loss_41(0)
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
        for (int i = 0; i < vertices_size; ++i)
        {
            loss_1 += vertices[i].state;
        }
        for (int i = 0; i < edge_size; ++i)
        {
            loss_21 += (vertices[edges[i][0]].state ^ vertices[edges[i][1]].state);
            loss_31 += (vertices[edges[i][0]].state ^ vertices[edges[i][1]].state) * (1.0L / vertices[edges[i][0]].size + 1.0L / vertices[edges[i][1]].size);
            loss_41 += (vertices[edges[i][0]].state ^ vertices[edges[i][1]].state) * (1.0L / (vertices[edges[i][0]].size * vertices[edges[i][0]].size) + 1.0L / (vertices[edges[i][1]].size * vertices[edges[i][1]].size));
        }
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
                if (w != u && !edge_matrix[u][w])
                {
                    // loss_21 += 2 * (vertices[v].state ^ vertices[w].state) - 1;
                    // loss_31 += (2 * (vertices[v].state ^ vertices[w].state) - 1) * (1.0 / vertices[v].size + 1.0 / vertices[w].size);
                    loss_41 += (2 * (vertices[v].state ^ vertices[w].state) - 1) * (1.0L / (vertices[v].size * vertices[v].size) + 1.0L / (vertices[w].size * vertices[w].size));
                }
            }
        }
    }

    long double loss(int choose)
    {
        if (choose == 1)
            return loss_1;
        else if (choose == 2)
            return loss_1 + temp_p * loss_21;
        else if (choose == 3)
            return loss_1 + temp_p * loss_31;
        else if (choose == 4)
            return loss_1 + temp_p * loss_41;
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
} g(0);

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
    // p_momentum = 0.9 * p_momentum + 10 * (temp_loss - loss) / (loss * T);
    loss = temp_loss;
    // T *= 1 - 1e-4;
    // num = 0;
    // temp_p += p_momentum;
    // if (temp_p > p_upbound)
    // {
    //     temp_p = p_upbound;
    //     p_momentum = 0;
    // }
    // if (temp_p < p_lowbound)
    // {
    //     temp_p = p_lowbound;
    //     p_momentum = 0;
    // }
}

void main_algorithm()
{
    g = graph(len);
    auto start_time = chrono::high_resolution_clock::now();
    loss = 100000;
    // loss = g.loss(loss_type);
    num = 0;
    sum_num = 0;
    temp_p = p;
    T = T_0;
    while (true)
    {
        // loss = g.loss(loss_type);
        int x = dis(gen) % len;
        g.light_up(x);
        temp_loss = g.loss(loss_type);
        if (loss >= temp_loss)
        {
            // T *= 1 + 0.2 * (temp_loss - loss) / loss;
            update();
        }
        else if (exp((loss - temp_loss) / T) > (dis_real(gen)))
        {
            // T *= 1 + 0.1 * (temp_loss - loss) / loss;
            update();
        }
        else
            g.light_up(x);
        num++;
        sum_num++;
        T *= cooling_rate;
        if (T < T_min)
        {
            break;
        }
        if (T < 1 && T > 0.1)
        {
            temp_p *= 1 - 2e-5;
        }
        if (T < 0.1)
        {
            temp_p *= 1 - 2e-5;
        }
        if (sum_num % 10000 == 0)
            cout << "temp_p: " << temp_p << "\tloss_num: " << g.loss_1 << "\tloss: " << loss << endl;
    }
    loss_num = g.loss_1;
    gradient_descent();
    loss = g.loss(loss_type);
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> runtime = end_time - start_time;
    // g.print();
    cout << "len: " << len << "\tp: " << p << "\tloss_num: " << loss_num << "\tloss: " << loss << " \tsum_num: " << sum_num << " \truntime: " << runtime.count() << endl;
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
    extract_last_entry(init_len, init_p);

    ofstream outfile("light_output.txt", std::ios::app);
    if (!outfile)
    {
        cerr << "无法打开文件 light_output.txt" << endl;
        exit(1);
    }

    for (len = len_min; len <= len_max; len += len_step)
    {
        for (p = p_min; p <= p_max + 1e-6; p += p_step)
        {
            if (is_init)
            {
                p = init_p + p_step;
                is_init = false;
                if (p > p_max)
                    break;
            }
            repeat_algorithm();
            outfile << "len: " << len << "\tp: " << p << "\taverage loss_num: " << sum / repeat_times << "\taverage runtime: " << sum_time / repeat_times;
            if (is_stable)
            {
                outfile << endl;
            }
            else
            {
                outfile << "\tnot stable: " << (sqrt(sum2 / (repeat_times * (repeat_times - 1)) - ((sum * sum) / (repeat_times * repeat_times * (repeat_times - 1)))) / (sum / repeat_times)) << endl;
            }
        }
    }

    outfile.close();
}

int main()
{
    // main_function();

    graph g(len);
    for (int i = 0; i < 100; ++i)
    {
        
    }

    return 0;
}
