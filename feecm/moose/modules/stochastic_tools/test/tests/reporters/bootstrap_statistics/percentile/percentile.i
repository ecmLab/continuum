[StochasticTools]
[]

[VectorPostprocessors]
  # Table 2.1 in Efron and Tibshirani, 1993
  [treatment]
    type = ConstantVectorPostprocessor
    value = '94 197 16 38 99 141 23'
    outputs = none
  []
  [control]
    type = ConstantVectorPostprocessor
    value = '52 104 146 10 51 30 40 27 46'
    outputs = none
  []
[]

[Reporters]
  # Reproduce Table 13.1 in Efron and Tibshirani, 1993
  [stats]
    type = StatisticsReporter
    vectorpostprocessors = 'treatment control'
    compute = 'mean stderr'
    ci_method = 'percentile'
    ci_levels = '0.025 0.05 0.1 0.16 0.5 0.84 0.9 0.95 0.975'
    ci_replicates = 1000
    ci_seed = 1980
  []
[]

[Outputs]
  execute_on = FINAL
  [out]
    type = JSON
  []
[]
