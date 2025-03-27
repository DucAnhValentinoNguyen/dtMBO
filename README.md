# dtMBO
Customed acquisition functions applying decision theory

Main files:
+ InfillFunctions.R: contains all the alternative functions
+ BenchmarkTess.R: implements the benchmark test

The test functions: 4 function classes: Ackley, Deflected Corrugated Srping, Rosenbrock, Schwefel, each has 3 functions (in 3 dimensions: 2D, 4D, 6D)
## Benchmark Test Functions

| Test Function | Formula | Properties | Domain | Minimum Location |
|--------------|---------|------------|--------|------------------|
| **Ackley** | `f(x) = −20exp(−0.2√(1/d)∑(x_i)) − exp((1/d)cos(2πx_i))` | Differentiable, multimodal | [-32.8, 32.8] | x<sub>i</sub> = 0 and f<sub>min</sub> = 0 |
| **Deflected Corrugated Spring** | `f(x) = 0.1∑(x_i - α)² - cos(K√∑(x_i - α)²)`<br>α = K = 5 by default | Highly multimodal, symmetric, global optimum far from local optima | [0, 2α] | x<sub>i</sub> = α and f<sub>min</sub> = −1 |
| **Rosenbrock** | `f(x) = ∑[100(x_{i+1} - x_i²)² + (1 - x_i)²]` | Unimodal, non-convex, differentiable | [-30, 30] | x<sub>i</sub> = 1 and f<sub>min</sub> = 0 |
| **Schwefel** | `f(x) = ∑[-x_i sin(√\|x_i\|)]` | Highly multimodal, global optimum far from local optima | [-500, 500] | x<sub>i</sub> = 420.9687<br>d=2: f<sub>min</sub> = −418.9829<br>d=4: f<sub>min</sub> = −1675.9316<br>d=6: f<sub>min</sub> = −2513.8974 |

Results:
![Benchmark Results](https://raw.githubusercontent.com/DucAnhValentinoNguyen/dtMBO/main/results/BenchmarkOutcome.png)
