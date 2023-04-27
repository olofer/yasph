#
# Read two diagnostics trace CSV files produced by YASPH.
# Require that the numbers in the same-size, same-column-name tables are within atol and rtol
#

import numpy as np 
import argparse 

def load_csv_and_colnames(csvname):
  with open(csvname) as f:
    first_line = f.readline().strip()
  Xnames = [e.strip() for e in first_line[1:].split(',')]
  X = np.loadtxt(csvname, skiprows = 1, delimiter = ',')
  assert X.shape[1] == len(Xnames)
  return X, Xnames

if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('--trace-1', type = str, default = '', help = 'name of CSV trace file 1')
  parser.add_argument('--trace-2', type = str, default = '', help = 'name of CSV trace file 2')

  parser.add_argument('--rtol', type = float, default = 1.0e-8, help = 'relative tolerance')
  parser.add_argument('--atol', type = float, default = 1.0e-4, help = 'relative tolerance')

  args = parser.parse_args()

  assert len(args.trace_1) > 0
  assert len(args.trace_2) > 0

  X1, Xnames1 = load_csv_and_colnames(args.trace_1)
  X2, Xnames2 = load_csv_and_colnames(args.trace_2)

  assert np.all(X1.shape == X2.shape)
  assert np.all([Xnames1[i] == Xnames2[i] for i in range(len(Xnames1))])

  maxabserr = np.max(np.abs(X1 - X2), axis = 0)
  maxabs1 = np.max(np.abs(X1), axis = 0)

  print(Xnames1)
  print(maxabserr)
  relerr = maxabserr / (1.0 + maxabs1)
  print(relerr)

  assert np.max(relerr) < args.rtol
  assert np.max(maxabserr) < args.atol 
