This is a pure vanilla implementation of LARS algorithm with LASSO modification.
For more information about the algorithm follow the link.

http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf

This document is included in the zip file (LeastAngle_2002.pdf).

There are many 'lars' implementation in the internet if you search. This is a relatively pure porting of the 'lars' package written by algorithm developer long time ago for R statistics language. Recent version by authors may have many changes. See http://cran.r-project.org/web/packages/lars/index.html for more information.

I commented equation number line by line in the codes so that you can follow the algorithm as you read the paper. You may start developing your own 'lars' implementation from this point. Test set includes Diabetes example used in the original paper. If you find any error, let me know. I'll appreciate it.

Caution: I've been using this file for more than 3 years, but I've never used 'forward stagewise' or 'pure lars'. I've been using only 'lars with lasso modification'. So I don't guarantee former two algorithms.


Installation:
Just copy 4 files in the 'lars' folder to your path. Or add 'lars' folder to your Matlab path.

Testing:
Read the 'test_lars.m' file. All files in the 'lars_tests' folder should be in the search path. It is written in the XTargets' MUnit format. But if you don't use MUnit, that's fine. Just copy and run codes in nested functions in that file.

Usage:
See 'test_lars.m' file.


License: The original R-package written by algorithm developers is under GPL v2 or newer. So is this.


