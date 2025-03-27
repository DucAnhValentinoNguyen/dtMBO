# dtMBO
Customed acquisition functions applying decision theory

Main files:
+ InfillFunctions.R: contains all the alternative functions
+ BenchmarkTess.R: implements the benchmark test

The test functions: 4 function classes: Ackley, Deflected Corrugated Srping, Rosenbrock, Schwefel, each has 3 functions (in 3 dimensions: 2D, 4D, 6D)
## Benchmark Test Functions

| Test Function       | Formula | Properties | Domain | Minimum Location |
|---------------------|---------|------------|--------|------------------|
| **Ackley** | \( f(x) = -20 \exp\left(-0.2 \sqrt{\frac{1}{d}\sum_{i=1}^d x_i^2}\right) - \exp\left(\frac{1}{d} \sum_{i=1}^d \cos(2\pi x_i)\right) + 20 + e \) | Differentiable, multimodal | \([-32.8, 32.8]\) | \(x_i = 0\), \(f_{\text{min}} = 0\) |
| **Deflected Corrugated Spring** | \( f(x) = 0.1 \sum_{i=1}^d (x_i - \alpha)^2 - \cos\left(K \sqrt{\sum_{i=1}^d (x_i - \alpha)^2}\right) \) <br> \(\alpha = K = 5\) (default) | Highly multimodal, symmetric, global optimum far from local optima | \([0, 2\alpha]\) | \(x_i = \alpha\), \(f_{\text{min}} = -1\) |
| **Rosenbrock** | \( f(x) = \sum_{i=1}^{d-1} \left[100(x_{i+1} - x_i^2)^2 + (1 - x_i)^2\right] \) | Unimodal, non-convex, differentiable | \([-30, 30]\) | \(x_i = 1\), \(f_{\text{min}} = 0\) |
| **Schwefel** | \( f(x) = \sum_{i=1}^d -x_i \sin\left(\sqrt{\|x_i\|}\right) \) | Highly multimodal, global optimum far from local optima | \([-500, 500]\) | \(x_i = 420.9687\) <br> \(d=2: f_{\text{min}} = -418.9829\) <br> \(d=4: f_{\text{min}} = -1675.9316\) <br> \(d=6: f_{\text{min}} = -2513.8974\) |

Results:
![Benchmark Results](https://raw.githubusercontent.com/DucAnhValentinoNguyen/dtMBO/main/results/BenchmarkOutcome.png)
