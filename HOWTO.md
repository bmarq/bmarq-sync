Evaluation of the bmarq-sync sycnhronization mechanism

usage: #./bmarq-sync-eval.py	[-h] [--version] [--nsim NSIM]
                         				[--maxcycles MAXCYCLES]
                        				[--alpha {0.01,0.05,0.1,0.125,0.25,0.375,0.5,0.625,0.75,0.875,0.9,0.95,0.99}]
                          				[--beta {1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10}]
                          				[--gamma {0.75,1.0,0.6,0.55,0.8,0.65,0.85,0.4,0.9,0.7,0.95,0.45,0.5}]
                          				[--sigma {0.25,0.0,0.1,0.05,0.15,0.2,0.01}]
                          				[--ton TON] [--tmaxcycle TMAXCYCLE]
                          				[--tmincycle TMINCYCLE] [--discard DISCARD]
                          				[--delaydist {constant,poisson,exponential,normal,uniform,chisquare}]
                          				[--cycledist {constant,poisson,exponential,normal,uniform,chisquare}]
                          				[--rndseed {True,False,T,F}]
                          
Evaluation of the bmarq-sync nodes sycnhronization mechanism

optional arguments:

  -h, --help			show this help message and exit

  --version             	Version number of eval_bmarq.py

  --nsim NSIM           	Number of simulations to perform (default: 1)

  --maxcycles 		MAXCYCLES
                        		Maximum number of cycles per simulation (default: 100)

  --alpha 			{0.01,0.05,0.1,0.125,0.25,0.375,0.5,0.625,0.75,0.875,0.9,0.95,0.99}
                        		Value for alpha parameter (default: 0.125)

  --beta 			{1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10}
                        		Value for beta parameter (default: 10)

  --gamma 		{0.75,1.0,0.6,0.55,0.8,0.65,0.85,0.4,0.9,0.7,0.95,0.45,0.5}
                        		% of TON for TSensorsOn success (default: 0.80)

  --sigma 			{0.25,0.0,0.1,0.05,0.15,0.2,0.01}
                        		Value of standard deviation for delays (default: 0.20)

  --ton TON             	Value for TON (default: 45)

  --tmaxcycle 		TMAXCYCLE
                        		Maximum value for TCycle (default: 3600)

  --tmincycle 		TMINCYCLE
                        		Minimum value for TCycle (default: 60)

  --discard 		DISCARD     Initial % of cycles to discard for estability purposes
                        		[default: 0.10 (10%)]

  --delaydist 		{constant,poisson,exponential,normal,uniform,chisquare}
                        		Type of random distribution for delays (default: uniform)

  --cycledist 		{constant,poisson,exponential,normal,uniform,chisquare}
                        		Type of random distribution for TCycle (default: uniform).
				If the distribution is constant, the default value equals TMINCYCLE

  --rndseed 		{True,False,T,F}
                        		Use prededined random seeds to reproduce experiments
                        		(T)rue/(F)alse (default: (F)alse)
