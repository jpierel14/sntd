from pylab import plot, show
fil = open('out_line')
all = fil.readlines()
x = []
y = []
for l in all:
    s = l.split()
    x.append(int(s[0]))
    y.append(float(s[4]))

plot(x,y)
show()


