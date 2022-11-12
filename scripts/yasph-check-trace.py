#
# Read diagnostics trace CSV file produced by YASPH.
# Check selected columns for conservation properties.
# Optionally write PDF figure files for each trace.
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

  parser.add_argument('--trace-file', type = str, default = '', help = 'name of CSV trace file')
  parser.add_argument('--trace-names', type = str, nargs = '+', help = 'name of columns in CSV file to be analyzed')
  parser.add_argument('--tolerance', type = float, default = None, help = '(optional) requirement for stdev / mean')
  parser.add_argument('--force-crash', action = 'store_true')
  parser.add_argument('--make-pdfs', action = 'store_true')

  args = parser.parse_args()

  assert len(args.trace_file) > 0

  X, Xnames = load_csv_and_colnames(args.trace_file)

  print(X.shape)
  print(Xnames)

  if args.make_pdfs:
    import matplotlib.pyplot as plt

  if not args.trace_names is None:
    print(args.trace_names)
    for xname in args.trace_names:
      assert xname in Xnames
      idx = Xnames.index(xname)
      x = X[:, idx]
      mean_x = np.mean(x)
      stdv_x = np.std(x)
      print('<{}> = {}, sd = {}'.format(xname, mean_x, stdv_x))

      if args.make_pdfs:
        plt.plot(X[:, Xnames.index('time')], x)
        plt.grid(True)
        plt.xlabel('time [sec]')
        plt.ylabel(xname)
        plt.title('<{}> = {:.4e}, sd = {:.4e}'.format(xname, mean_x, stdv_x))
        plt.savefig('yasph-check-trace-{}.pdf'.format(xname))
        plt.close()

      if args.tolerance is not None:
        ratio = stdv_x / np.abs(mean_x) if np.abs(mean_x) > 1 else stdv_x
        print('ratio = {} (ok={})'.format(ratio, ratio < args.tolerance))
        if args.force_crash:
          assert ratio < args.tolerance

  else:
    print('no trace names specified')
