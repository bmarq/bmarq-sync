Evaluation of the bmarq-sync sycnhronization mechanism

usage: eval_bmarq-sync.py [-h] [--version] [--nsim NSIM]
                          [--maxcicles MAXCICLES] [--alpha {0.125,0.5,0.875}]
                          [--beta {1,10,50,100}]
                          [--gamma {0.5,1.0,0.8,0.6,0.85,0.9,0.95,0.7}]
                          [--sigma {0.25,0.0,0.1,0.05,0.15,0.2}] [--TON TON]
                          [--tMaxCycle TMAXCYCLE] [--tMinCycle TMINCYCLE]
                          [--pDiscard PDISCARD]
                          [--tDelayDist {normal,uniform,chisquare,poisson,exponential}]
                          [--tCycleDist {normal,uniform,chisquare,poisson,exponential}]
                          
optional arguments:
  -h, --help            	show this help message and exit
  --version             	Version number of eval_bmarq.py
  --nsim NSIM           	Number of simulations to perform (default: 1)
  --maxcicles MAXCICLES		Maximum number of cycles per simulation (default: 1000)
  --alpha {0.125,0.5,0.875}	Value for alpha parameter default: 0.125)
  --beta {1,10,50,100}		Value for beta parameter (default: 10)
  --gamma {0.5,1.0,0.8,0.6,0.85,0.9,0.95,0.7}	% of TON for TSensorsOn success (default: 0.80)
  --sigma {0.25,0.0,0.1,0.05,0.15,0.2}	Value of standard deviation for delays (default: 0.20)
  --TON TON             	Value for TON (default: 60)
  --tMaxCycle TMAXCYCLE		Maximum value for TCycle (default: 900)
  --tMinCycle TMINCYCLE		Minimum value for TCycle (default: 120)
  --pDiscard PDISCARD   Initial % of cycles to discard for estability purposes	[default: 0.10 (10%)]
  --tDelayDist {normal,uniform,chisquare,poisson,exponential}	Type of random distribution for delays (default: uniform)
  --tCycleDist {normal,uniform,chisquare,poisson,exponential}	Type of random distribution for TCycle (default: uniform)
