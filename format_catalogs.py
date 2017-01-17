import numpy
import readfile

data = readfile.table('udgs_vdburg16.dat', dtype=(str,int,float,float,float))
clusters = set(data[0])

for cluster in clusters:
    j = numpy.arange(data[0].size)[data[0] == cluster]
    print cluster, len(j)
    catalog = 'catalogs/{0}.cat'.format(cluster)
    cat = open(catalog, 'w')
    print >>cat, '# id  x  y  rmag'
    for i in j:
        obj = [d[i] for d in data[1:]]
        print >>cat, '{0:6d}  {1:9.3f}  {2:9.3f}  {3:5.2f}'.format(*obj)
    cat.close()
