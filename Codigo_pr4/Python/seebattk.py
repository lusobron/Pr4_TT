def seebattk(v):
    # Locals
    d = [1.0/3.0, 4.0/27.0, 8.0/27.0, 2.0/9.0, 22.0/81.0, 208.0/891.0, 
         340.0/1287.0, 418.0/1755.0, 598.0/2295.0, 700.0/2907.0, 928.0/3591.0, 
         1054.0/4347.0, 1330.0/5175.0, 1480.0/6075.0, 1804.0/7047.0, 
         1978.0/8091.0, 2350.0/9207.0, 2548.0/10395.0, 2968.0/11655.0, 
         3190.0/12987.0, 3658.0/14391.0]
    
    # Process Forwards
    sum1 = d[0]
    delold = 1.0
    termold = d[0]
    i = 0
    while True:
        del_ = 1.0 / (1.0 + d[i+1] * v * delold)
        term = termold * (del_ - 1.0)
        sum1 = sum1 + term
        i = i + 1
        delold = del_
        termold = term
        if i > 20 or abs(termold) > 0.000001:
            break
    
    return sum1