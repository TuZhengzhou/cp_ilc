#include <iostream>
#include <cmath>
#include <limits>

// Objective function
double objectiveFunction(int n1_, int n2_, int m_) {
    return (n1_ + n2_) * log(n1_ + n2_) + (m_ + m_) * (n1_ + n2_);
}

unsigned int nextPowerOf2(unsigned int n) {
    // 如果 n 已经是2的幂次方，则直接返回
    if ((n & (n - 1)) == 0) {
        return n;
    }

    unsigned int count = 0;

    // 如果 n 不是2的幂次方，则计算向上取整的2的幂次方值
    while (n != 0) {
        n >>= 1;
        count++;
    }

    return 1 << count;
}

// Function to check if a number is a power of 2
bool isPowerOfTwo(int num) {
    return (num != 0) && ((num & (num - 1)) == 0);
}

void solveOptimizationProblem(double t1, double t2, size_t& n1, size_t& n2, size_t& m) {
    double min_result = std::numeric_limits<double>::max();
    int opt_n1 = -1, opt_n2 = -1, opt_m = -1;
    int n1_, n2_, m_;
    unsigned int m_upper;
    double t1_div_n1, t2_div_n2;
    for (n1_ = 1; n1_ < t1; ++n1_) {
        for (n2_ = 1; n2_ < t2; ++n2_) {
            t1_div_n1 = t1 / n1_;
            t2_div_n2 = t2 / n2_;
            m_upper = nextPowerOf2(std::ceil(std::max(t1_div_n1, t2_div_n2)));
            for (m_ = 1; m_ <= m_upper; ++m_) {
                if (m_ * n1_ >= t1 && m_ * n2_ >= t2 && isPowerOfTwo(m_)) {
                    double result = objectiveFunction(n1_, n2_, m_);
                    if (result < min_result) {
                        min_result = result;
                        opt_n1 = n1_;
                        opt_n2 = n2_;
                        opt_m = m_;
                    }
                }
            }
        }
    }
    
    n1 = opt_n1;
    n2 = opt_n2;
    m  = opt_m;
    std::cout << "Optimal values for n1_, n2_, m_: " << opt_n1 << ", " << opt_n2 << ", " << opt_m << std::endl;
}

int main() {
    int t1 = 300, t2 = 350; // Given inputs t1, t2
    size_t n1, n2, m;
    solveOptimizationProblem(t1, t2, n1, n2, m);

    return 0;
}
