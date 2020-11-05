# Learn-to-Compress
To have a quick start, you can run to generate data needed
```
python gen_norm.py

```
Then, you can run the command
```
g++ piecewise.cpp -o piecewise
./piecewise -f filename -e max_error
```
In which there are two hyper-parameter, filename in {linear,noisylinear,normal,noisynormal,lognormal,books}
and max_error is an arbitrary int you choose (2^k - 1 like 7,15,31,63,127... is recommended)
