import time

def dtimer(tim, itime, icntrl):
# input: icntrl, (float) itime
# icntrl = (-1,0,1) = (initialize,ignore,read) clock
# clock should be initialized before it is read!
# returns (float) tim[0] = elapsed time in seconds
    if icntrl == 0:
        return
    if icntrl == 1:
        usec = time.time() - itime[0]
        tim[0] = usec
    else:
        itime[0] = time.time()
