# SafePatternPruning: An Efficient Approach for Predictive Pattern Mining (KDD'16)

## Usage

`./train [option] [filename]`

### option: 
- -T : compute regularization path for a sequence of T \lambda evenly allocated between \lambda_0 and 0.01\lambda

- -F : calculate duality gap and dynamic screening every F iteration

- -S : minimum support (graph only)

- -D : maxdepth or maxpat

- -B : with bias (0 or 1)
