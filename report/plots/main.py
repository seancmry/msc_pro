import matplotlib.pyplot as plt

plt.style.use('seaborn-whitegrid')

sd_ackley = [ 0.11635, 0.17158, 0.29957, 0.28402, 0.41190 ]
sd_sphere = [ 40.254, 146.13, 146.14, 243.31, 435.64 ]
sd_rosenbrock = [ 155.74, 322.99, 1942.9, 1951.6, 3699.4 ]
sd_griewank = [ 1.022, 0.76188, 1.1489, 1.2648, 1.3156 ]
dims = [ 20, 30, 50, 70, 100 ]

sd_ackley2 = [ 0.17158, 0.41190, 0.61797, 0.32277 ]
sd_sphere2 = [ 146.13, 435.64, 27717.0, 112470.0 ]
sd_rosenbrock2 = [ 322.99, 3699.4, 6974.8, 10289.0 ]
sd_griewank2 = [ 0.76188, 1.3156, 71.154, 328.38 ]
dims2 = [ 30, 100, 500, 1000 ]

step = [ 0,
         50,
         100,
         150,
         200,
         250,
         300,
         350,
         400,
         450,
         500,
         550,
         600,
         650,
         700,
         750,
         800,
         850,
         900,
         950,
         1000,
         1050,
         1100,
         1150,
         1200,
         1250,
         1300,
         1350,
         1400,
         1450,
         1500,
         1550,
         1600,
         1650,
         1700,
         1750,
         1800,
         1850,
         1900,
         1950,
         2000,
         2050,
         2100,
         2150,
         2200,
         2250,
         2300,
         2350,
         2400,
         2450,
         2500,
         2550,
         2600,
         2650,
         2700,
         2750,
         2800,
         2850,
         2900,
         2950,
         3000,
         3050,
         3100,
         3150,
         3200,
         3250,
         3300,
         3350,
         3400,
         3450,
         3500,
         3550,
         3600,
         3650,
         3700,
         3750,
         3800,
         3850,
         3900,
         3950,
         4000,
         4050,
         4100,
         4150,
         4200,
         4250,
         4300,
         4350,
         4400,
         4450,
         4500,
         4550,
         4600,
         4650,
         4700,
         4750,
         4800,
         4850,
         4900,
         4950,
         5000,
         5050,
         5100,
         5150,
         5200,
         5250,
         5300,
         5350,
         5400,
         5450,
         5500,
         5550,
         5600,
         5650,
         5700,
         5750,
         5800,
         5850,
         5900,
         5950,
         6000,
         6050,
         6100,
         6150,
         6200,
         6250,
         6300,
         6350,
         6400,
         6450,
         6500,
         6550,
         6600,
         6650,
         6700,
         6750,
         6800,
         6850,
         6900,
         6950,
         7000,
         7050,
         7100,
         7150,
         7200,
         7250,
         7300,
         7350,
         7400,
         7450,
         7500,
         7550,
         7600,
         7650,
         7700,
         7750,
         7800,
         7850,
         7900,
         7950,
         8000,
         8050,
         8100,
         8150,
         8200,
         8250,
         8300,
         8350,
         8400,
         8450,
         8500,
         8550,
         8600,
         8650,
         8700,
         8750,
         8800,
         8850,
         8900,
         8950,
         9000,
         9050,
         9100,
         9150,
         9200,
         9250,
         9300,
         9350,
         9400,
         9450,
         9500,
         9550,
         9600,
         9650,
         9700,
         9750,
         9800,
         9850,
         9900,
         9950,
         10000
         ]
mean_ackley = [ 1.98E+01,
                1.91E+01,
                1.88E+01,
                1.86E+01,
                1.83E+01,
                1.80E+01,
                1.77E+01,
                1.75E+01,
                1.72E+01,
                1.69E+01,
                1.67E+01,
                1.64E+01,
                1.62E+01,
                1.60E+01,
                1.57E+01,
                1.55E+01,
                1.53E+01,
                1.51E+01,
                1.48E+01,
                1.47E+01,
                1.45E+01,
                1.44E+01,
                1.42E+01,
                1.39E+01,
                1.37E+01,
                1.35E+01,
                1.34E+01,
                1.32E+01,
                1.30E+01,
                1.28E+01,
                1.27E+01,
                1.25E+01,
                1.23E+01,
                1.21E+01,
                1.19E+01,
                1.17E+01,
                1.15E+01,
                1.14E+01,
                1.13E+01,
                1.11E+01,
                1.09E+01,
                1.07E+01,
                1.06E+01,
                1.05E+01,
                1.04E+01,
                1.03E+01,
                1.01E+01,
                9.99E+00,
                9.85E+00,
                9.76E+00,
                9.70E+00,
                9.60E+00,
                9.50E+00,
                9.42E+00,
                9.35E+00,
                9.26E+00,
                9.19E+00,
                9.12E+00,
                9.02E+00,
                8.95E+00,
                8.88E+00,
                8.81E+00,
                8.74E+00,
                8.68E+00,
                8.64E+00,
                8.57E+00,
                8.51E+00,
                8.46E+00,
                8.41E+00,
                8.37E+00,
                8.32E+00,
                8.28E+00,
                8.26E+00,
                8.22E+00,
                8.18E+00,
                8.15E+00,
                8.12E+00,
                8.09E+00,
                8.05E+00,
                8.02E+00,
                8.00E+00,
                7.98E+00,
                7.97E+00,
                7.95E+00,
                7.93E+00,
                7.92E+00,
                7.91E+00,
                7.90E+00,
                7.90E+00,
                7.89E+00,
                7.88E+00,
                7.88E+00,
                7.88E+00,
                7.88E+00,
                7.88E+00,
                7.88E+00,
                7.88E+00,
                7.88E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.87E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00,
                7.86E+00
                ]
mean_sphere = [ 2.83E+05,
                1.47E+05,
                1.08E+05,
                8.26E+04,
                6.54E+04,
                5.42E+04,
                4.45E+04,
                3.66E+04,
                3.05E+04,
                2.56E+04,
                2.22E+04,
                1.87E+04,
                1.62E+04,
                1.45E+04,
                1.30E+04,
                1.17E+04,
                1.04E+04,
                9.44E+03,
                8.67E+03,
                7.89E+03,
                7.17E+03,
                6.59E+03,
                5.98E+03,
                5.52E+03,
                5.06E+03,
                4.69E+03,
                4.42E+03,
                4.12E+03,
                3.80E+03,
                3.57E+03,
                3.34E+03,
                3.11E+03,
                2.93E+03,
                2.76E+03,
                2.59E+03,
                2.43E+03,
                2.31E+03,
                2.17E+03,
                2.07E+03,
                1.96E+03,
                1.86E+03,
                1.80E+03,
                1.70E+03,
                1.63E+03,
                1.56E+03,
                1.50E+03,
                1.45E+03,
                1.39E+03,
                1.35E+03,
                1.30E+03,
                1.26E+03,
                1.21E+03,
                1.18E+03,
                1.15E+03,
                1.12E+03,
                1.09E+03,
                1.07E+03,
                1.04E+03,
                1.01E+03,
                9.89E+02,
                9.75E+02,
                9.56E+02,
                9.32E+02,
                9.20E+02,
                9.01E+02,
                8.85E+02,
                8.73E+02,
                8.57E+02,
                8.47E+02,
                8.32E+02,
                8.23E+02,
                8.12E+02,
                8.03E+02,
                7.94E+02,
                7.87E+02,
                7.80E+02,
                7.73E+02,
                7.66E+02,
                7.61E+02,
                7.56E+02,
                7.50E+02,
                7.45E+02,
                7.41E+02,
                7.36E+02,
                7.35E+02,
                7.32E+02,
                7.30E+02,
                7.28E+02,
                7.26E+02,
                7.25E+02,
                7.24E+02,
                7.23E+02,
                7.22E+02,
                7.21E+02,
                7.20E+02,
                7.20E+02,
                7.20E+02,
                7.19E+02,
                7.19E+02,
                7.19E+02,
                7.18E+02,
                7.18E+02,
                7.18E+02,
                7.18E+02,
                7.18E+02,
                7.17E+02,
                7.17E+02,
                7.17E+02,
                7.17E+02,
                7.17E+02,
                7.17E+02,
                7.17E+02,
                7.17E+02,
                7.16E+02,
                7.16E+02,
                7.16E+02,
                7.16E+02,
                7.16E+02,
                7.16E+02,
                7.16E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.15E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02,
                7.14E+02
                ]
mean_rosenbrock = [ 3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04,
                    3.50E+04
                    ]
mean_griewank = [ 2.54E+03,
                  1.16E+03,
                  7.75E+02,
                  5.10E+02,
                  3.62E+02,
                  2.63E+02,
                  1.92E+02,
                  1.40E+02,
                  1.02E+02,
                  7.71E+01,
                  5.85E+01,
                  4.55E+01,
                  3.47E+01,
                  2.81E+01,
                  2.17E+01,
                  1.70E+01,
                  1.38E+01,
                  1.13E+01,
                  9.21E+00,
                  7.84E+00,
                  6.73E+00,
                  5.78E+00,
                  5.08E+00,
                  4.50E+00,
                  3.99E+00,
                  3.67E+00,
                  3.36E+00,
                  3.08E+00,
                  2.83E+00,
                  2.65E+00,
                  2.50E+00,
                  2.37E+00,
                  2.24E+00,
                  2.14E+00,
                  2.05E+00,
                  1.98E+00,
                  1.91E+00,
                  1.84E+00,
                  1.79E+00,
                  1.74E+00,
                  1.69E+00,
                  1.64E+00,
                  1.61E+00,
                  1.57E+00,
                  1.54E+00,
                  1.51E+00,
                  1.48E+00,
                  1.46E+00,
                  1.43E+00,
                  1.41E+00,
                  1.39E+00,
                  1.37E+00,
                  1.36E+00,
                  1.35E+00,
                  1.33E+00,
                  1.32E+00,
                  1.32E+00,
                  1.31E+00,
                  1.30E+00,
                  1.29E+00,
                  1.28E+00,
                  1.27E+00,
                  1.26E+00,
                  1.25E+00,
                  1.25E+00,
                  1.24E+00,
                  1.24E+00,
                  1.23E+00,
                  1.23E+00,
                  1.23E+00,
                  1.22E+00,
                  1.22E+00,
                  1.22E+00,
                  1.21E+00,
                  1.21E+00,
                  1.21E+00,
                  1.21E+00,
                  1.21E+00,
                  1.21E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.20E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00,
                  1.19E+00
                  ]

time_ackley = [ 270.4, 2908.0, 147818.00, 560022.8 ]
time_sphere = [ 160.0, 1819.6, 95558.4, 419196.8 ]
time_rosenbrock = [ 240.8, 2831.2, 74689.6, 299853.2 ]
time_griewank = [ 290.4, 3194.0, 134804.8, 590098.0 ]

plt.plot(dims, sd_ackley, "-s", label="Ackley")
plt.plot(dims, sd_sphere, "-o", label="Sphere")
plt.plot(dims, sd_rosenbrock, "-v", label="Rosenbrock")
plt.plot(dims, sd_griewank, "-h", label="Griewank")
plt.xlabel("Dimensions (N)")
plt.ylabel("Standard deviation of aggregate means from individual runs")
plt.title("Standard deviation of mean best solution over 25 runs across N dimensions")
plt.legend(loc="upper left")
plt.show()

plt.plot(dims2, sd_ackley2, "-s", label="Ackley")
plt.plot(dims2, sd_sphere2, "-o", label="Sphere")
plt.plot(dims2, sd_rosenbrock2, "-v", label="Rosenbrock")
plt.plot(dims2, sd_griewank2, "-h", label="Griewank")
plt.xlabel("Dimensions (N)")
plt.ylabel("Standard deviation of aggregate means from individual runs")
plt.title("Standard deviation of mean best solution over 25 runs across N dimensions")
plt.legend(loc="upper left")
plt.show()

plt.plot(dims2, time_ackley, "-s", label="Ackley")
plt.plot(dims2, time_sphere, "-o", label="Sphere")
plt.plot(dims2, time_rosenbrock, "-v", label="Rosenbrock")
plt.plot(dims2, time_griewank, "-h", label="Griewank")
plt.xlabel("Dimensions (N)")
plt.ylabel("Mean time in milliseconds (ms)")
plt.title("Average time per run across N dimensions")
plt.legend(loc="upper left")
plt.show()


plt.plot(step, mean_ackley)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.xlabel("Iterations")
plt.ylabel("Mean fitness")
plt.title("Ackley")
plt.show()

plt.plot(step, mean_sphere)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.xlabel("Iterations")
plt.ylabel("Mean fitness")
plt.title("Sphere")
plt.show()

plt.plot(step, mean_rosenbrock)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.xlabel("Iterations")
plt.ylabel("Mean fitness")
plt.title("Rosenbrock")
plt.show()

plt.plot(step, mean_griewank)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.xlabel("Iterations")
plt.ylabel("Mean fitness")
plt.title("Griewank")
plt.show()