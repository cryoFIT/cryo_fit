##!/usr/bin/env python
#!/usr/bin/python

# When an warning/error message of
#/Users/ftuser/bin/phenix-dev-2880/base/Python.framework/Versions/2.7/lib/python2.7/site-packages/matplotlib-2.0.0-py2.7-macosx-10.6-x86_64.egg/matplotlib/axes/_axes.py:545: UserWarning: No labelled objects found. Use label='...' kwarg on individual plots.
#  warnings.warn("No labelled objects found. "
# appears,
'''
sudo easy_install pip
pip install matplotlib
didn't help (9/15/2017, lanl ftuser)
'''

# when 'matplotlib is not installed in mac'
'''
if you are on OS X 10.7 or 10.8, using the Apple-installed Python, you
have easy_install built-in, but not pip. To fix that:
   sudo easy_install pip
   And now, you can do this:  sudo pip install matplotlib
<ref> http://stackoverflow.com/questions/13979496/python-import-matplotlib-pyplot-not-working
'''

# when 'no display name and no $DISPLAY environment variable error' when user tries to run this cc_vs_score.py on linux instead of mac
'''
When signing into the server to execute the code use this instead:
   ssh -X username@servername
   the -X will get rid of the no display name and no $DISPLAY environment
   variable error
<ref> http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
'''

import math, os, random, sys
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 12,
        }

try:
   print "provided cor.out: ", sys.argv[1]
except:
   print "cor.out is not provided as an argument, so quit the program"
   sys.exit(0)

print "draw step vs cc plot"

cc_array = []
step_array = []

with open(sys.argv[1]) as f:
   for line in f:
      splited = line.split()
      step = splited[1]
      step_array.append(float(step))
      cc = splited[4]
      cc_array.append(float(cc))

print "min(cc_array): ", min(cc_array)
print "max(cc_array): ", max(cc_array)

print "min(step_array): ", min(step_array)
print "max(step_array): ", max(step_array)

input_file = sys.argv[1]

#plt.plot(cc_array, score_array, linestyle='--', marker='o', label='each decoy')
#plt.plot(cc_array, step_array, 'bs', marker='o', label='decoys')
plt.scatter(step_array, cc_array, 6, color='red')

plt.xlim([min(step_array)-1, max(step_array)+1])
plt.ylim([min(cc_array)-0.001, max(cc_array)+0.001])

plt.text(min(step_array)-100, max(step_array)+100, 0.87, fontdict=font)

plt.xlabel('step')
plt.ylabel('cc')
plt.title(input_file)
plt.legend()
plt.show()
