import pylab as pl
import scipy.stats as sps
import numpy.random as nr
import numpy as np
import networkx as nw
import collections as col

ext='.png'

"""--------------------------------PRICE MODEL CLASS--------------------------------"""
class BAmodel(object):
    def __init__(self,m,pref):#constructor
        self.m=m        #edges added per time increment
        self.pref=pref        #network "knows" whether it preferentially attaches
        if (pref):
            self.G=nw.complete_graph(m+1)
            self.deg_nodelist=[]#degenerate node list (degeneracy=degree)
            self.deg_nodelist.extend(range(m+1)*m)
        else:
            self.G=nw.empty_graph(m)

    def gephit(self,filename):
        """output as a gephi-readable file"""
        nw.write_gexf(self.G,filename)

    def evolve(self,time):        #evolve for a given time
        for i in range(time):
            newly_attached=self.choose_nodes()
            new_node=[self.G.number_of_nodes()]*self.m
            self.G.add_edges_from( zip( new_node , newly_attached ))            #add new node, connect m edges to chosen nodes
            if (self.pref):#update the list of nodes repeated for degree if the system is preferentially attaching
                self.deg_nodelist.extend(newly_attached)
                self.deg_nodelist.extend(new_node)

    def choose_nodes(self):
        chosen_nodes=set()        #sets are by definition unique so that the same node can't be rechosen
        if (self.pref):#preferential attachment case
            while(len(chosen_nodes)<self.m):
                chosen_nodes.add(nr.choice(self.deg_nodelist))#choose from the deg list of nodes
        else:#random attachment case
            while(len(chosen_nodes)<self.m):
                chosen_nodes.add(nr.choice(self.G.nodes()))#choose from the list of nodes
        return chosen_nodes

    def walk_evolve(self,time,L):
        for i in range(time):
            newly_attached=self.random_walk(L)
            new_node=[self.G.number_of_nodes()]*self.m
            self.G.add_edges_from(zip(new_node , newly_attached ))

    def random_walk(self,L):
        nodes=set()
        while(len(nodes)<self.m):
            b=nr.choice(self.G.nodes())#choose a node at random
            for i in range(L):
                b=nr.choice(self.G.neighbors(b))#choose nodes from a list of neighbours of previous
            nodes.add(b)#add chosen node to list.
        return nodes


"""-----------------------------------FUNCTIONS-------------------------------------"""

def degree_freqdist(degree_list,base=None):#nw.degree(graph).values()=degreelist
    freqdist={}
    if (base!=None):
        for i in degree_list:
            if (i!=0):
                bin_number=int(pl.log(i)/pl.log(base))
            else:
                bin_number=0
            if(bin_number in freqdist):
                freqdist[bin_number]+=1
            else:
                freqdist[bin_number]=1
    else:
        for i in degree_list:
            if (i in freqdist):
                freqdist[i]+=1
            else:
                freqdist[i]=1
    return freqdist
#for probability, divide by N=m+1+time (PA); N=m+time (rnd), or could just query the graph object

def binvalues(bin_number,base):
    #returns binvalue, bin min, bin max for given bin number and base
    bin_max=float(base)
    bin_min=1.0
    for i in range(bin_number):
        bin_min*=base
        bin_max*=base
    bin_min=pl.floor(bin_min)
    bin_max=pl.floor(bin_max)
    if (bin_min!=bin_max):
        bin_min+=1.0
    return (pl.sqrt(bin_min*bin_max),bin_min,bin_max)

def f_to_p(freq_dist, N, base=None):#converts frequency to probability
    prob_dist={}
    if (base!=None):
        for b_number,freq in freq_dist.items():
            (degree,d_min,d_max)=binvalues(b_number,base)
            prob_dist[degree]=float(freq)/((d_max-d_min+1.0)*float(N))
    else:
        for degree,freq in freq_dist.items():
            prob_dist[degree]=float(freq)/float(N)
    return prob_dist
    


"""-----PHASE 1:PART 3, DEGREE DISTRIBUTION-----"""
def phase13():
    nets=[]
    results=[]
    m=1
    it=100
    av=100
    base=1.1
    f = open("chi-square.txt", "w")
    f.write('m\tchi_square\tp_value\td.o.f.\n')
    for i in range(4):
        res={}
        nets.append(BAmodel(m,True))
        for j in range(av):
            interes={}
            nets[i]=BAmodel(m,True)
            nets[i].evolve(it)
            interes=degree_freqdist( nw.degree(nets[i].G).values() , base )
            for degree, freq in interes.items():
                if (degree in res):
                    res[degree]+=freq
                else:
                    res[degree]=freq
            print '%dth of 4 Graphs %d %% done' % (i+1,(j+1))
        for degree in res.keys():
            res[degree]/=av
        results.append(col.OrderedDict(sorted( f_to_p(res, it+m+1, base).items() )))
        observed=pl.array(results[i].values())*(it+m+1)#observed frequency
        exp=pl.array(results[i].keys())
        exp=(2*m*(m+1)/(exp*(exp+1)*(exp+2)))#expected probability
        pl.loglog()
        pl.ylabel('p(k)')
        pl.xlabel('k')
        pl.title('Degree PDF for m='+str(m) )
        pl.plot(results[i].keys(),results[i].values(), 'r-')
        pl.plot(results[i].keys(), exp, 'b--')
        pl.legend(('Data', 'Theoretical Form'),loc='upper right')
        pl.savefig('comparison'+str(i)+ext)
        pl.clf()
        exp*=(it+m+1)#expected frequency
        (chisq,p_value)=sps.chisquare(observed,exp)
        f.write(str(m)+'\t'+str(chisq)+'\t'+str(p_value)+'\t'+str(len(exp)-1)+'\n')
        m*=2

    pl.loglog()
    pl.ylabel('p(k)')
    pl.xlabel('k')
    pl.title('Degree PDF' )
    pl.plot(results[0].keys(), results[0].values(), 'r-')
    pl.plot(results[1].keys(), results[1].values(), 'b-')
    pl.plot(results[2].keys(), results[2].values(), 'g-')
    pl.plot(results[3].keys(), results[3].values(), 'm-')
    pl.legend(('m=1', 'm=2', 'm=4', 'm=8'),loc='upper right')
    pl.savefig('degree_dist'+ext)
    pl.clf()

    f.close()

"""-----PHASE 1:PART 4, FINITE SIZE EFFECT------"""
def phase14():
    f = open("largest k.dat", "w")
    f1= open("kldata.dat","w")
    f.write('N\tlargest k\tpredicted largest k\terror\n')
    nets=[]
    results=[]
    largest_k=[]
    mslargestk=[]
    m=8
    it=100
    av=[1000,100,5]
    base=1.1
    pl.loglog()
    pl.ylabel('p(k)')
    pl.xlabel('k')
    pl.title('Degree PDF for different number of nodes' )
    for i in range(3):
        largest_k.append(0.0)
        mslargestk.append(0.0)
        res={}
        nets.append(BAmodel(m,True))
        for j in range(av[i]):
            interes={}
            nets[i]=BAmodel(m,True)
            nets[i].evolve(it)
            degrees=sorted( nw.degree(nets[i].G).values())
            largest_k[i]+=degrees[len(degrees)-1]
            mslargestk[i]+=degrees[len(degrees)-1]*degrees[len(degrees)-1]
            f1.write( str(i) +'\t'+str(degrees[len(degrees)-1])+'\t'+str(degrees[len(degrees)-1]*degrees[len(degrees)-1])+'\n')
            interes=degree_freqdist(degrees , base )
            for degree, freq in interes.items():
                if (degree in res):
                    res[degree]+=freq
                else:
                    res[degree]=freq
            print '%dth of 3 Graphs %d %% done' % ((i+1),float((j+1)*100)/float(av[i]))
        for degree in res.keys():
            res[degree]/=av[i]
        largest_k[i]/=av[i]
        mslargestk[i]/=av[i]
        f.write(str(it)+'\t'+str(largest_k[i])+'\t'+str(pl.sqrt(2*(it+m)))+'\t'+str(pl.sqrt((mslargestk[i]-largest_k[i]*largest_k[i]))/av[i])+'\n')
        results.append(col.OrderedDict(sorted( f_to_p(res, it+m, base).items() )))
        it*=10
    pl.plot(results[0].keys(), results[0].values(), 'r-')
    pl.plot(results[1].keys(), results[1].values(), 'b-')
    pl.plot(results[2].keys(), results[2].values(), 'g-')
    #pl.plot(results[3].keys(), results[3].values(), 'm-')
    pl.legend(('N=100', 'N=1000', 'N=10000', 'N=100000'),loc='upper right')
    pl.savefig('degree_dist_varying_N'+ext)
    pl.clf()

    pl.loglog()
    pl.ylabel('p(k)/(p_th(k))')
    pl.xlabel('k/k_l')
    pl.title('Data collapsed pdfs' )

    obs=pl.array(results[0].values())
    val=pl.array(results[0].keys())
    exp=(2*m*(m+1)/(val*(val+1)*(val+2)))
    pl.plot(val/largest_k[0], obs/exp, 'r-')

    obs=pl.array(results[1].values())
    val=pl.array(results[1].keys())
    exp=(2*m*(m+1)/(val*(val+1)*(val+2)))
    pl.plot(val/largest_k[1], obs/exp, 'b-')

    obs=pl.array(results[2].values())
    val=pl.array(results[2].keys())
    exp=(2*m*(m+1)/(val*(val+1)*(val+2)))
    pl.plot(val/largest_k[2], obs/exp, 'g-')

    """obs=pl.array(results[3].values())
    exp=pl.array(results[3].keys())
    exp=(2*m*(m+1)/(val*(val+1)*(val+2)))
    pl.plot(results[3].keys(), results[3].values(), 'm-')"""
    pl.legend(('N=100', 'N=1000', 'N=10000', 'N=100000'),loc='upper right')
    pl.savefig('degree_dist_collapsed'+ext)
    pl.clf()

    f.close()
    f1.close()


        
"""-----PHASE 2, DEGREE DISTRIBUTION:PRA--------"""

def phase2():
    nets=[]
    results=[]
    largest_k=[]
    m=1
    it=10000
    av=100
    f = open("chi-square_pra.txt", "w")
    f1= open("largest k.dat","w")
    f1.write('m\tlargest k\tpredicted largest k\n')
    f.write('m\tchi_square\tp_value\td.o.f.\n')
    for i in range(4):
        largest_k.append(0.0)
        res={}
        nets.append(BAmodel(m,False))
        for j in range(av):
            interes={}
            nets[i]=BAmodel(m,False)
            nets[i].evolve(it)
            values=sorted( nw.degree(nets[i].G).values())
            largest_k[i]+=float(values[len(values)-1])
            interes=degree_freqdist(values)
            for degree, freq in interes.items():
                if (degree in res):
                    res[degree]+=freq
                else:
                    res[degree]=freq
            print '%dth of 4 Graphs %d %% done' % (i+1,(j+1))
        for degree in res.keys():
            res[degree]/=av
        largest_k[i]/=av
        results.append(col.OrderedDict(sorted( f_to_p(res, it+m).items() )))
        observed=pl.array(results[i].values())*(it+m)#observed frequency
        exp=pl.array(results[i].keys())
        exp=pl.array(float(m)-exp)
        C=(float(m+1)/float(m))
        exp=np.power( C , exp  )#expected probability
        exp/=(m+1)
        pl.loglog()
        pl.ylabel('p(k)')
        pl.xlabel('k')
        pl.title('Degree PDF for m='+str(m) )
        pl.plot(results[i].keys(),results[i].values(), 'r-')
        pl.plot(results[i].keys(), exp, 'b--')
        pl.legend(('Data', 'Theoretical Form'),loc='upper right')
        pl.savefig('comparison'+str(i)+ext)
        pl.clf()
        exp*=(it+m)#expected frequency
        (chisq,p_value)=sps.chisquare(observed,exp)
        f.write(str(m)+'\t'+str(chisq)+'\t'+str(p_value)+'\t'+str(len(exp)-1)+'\n')
        A=pl.log(C)
        k_largest=(1/A)*pl.log((C**(float(m)))*(it+m))
        f1.write(str(m)+'\t'+str(largest_k[i])+'\t'+str(k_largest)+'\n')
        m*=2

    pl.loglog()
    pl.ylabel('p(k)')
    pl.xlabel('k')
    pl.title('Degree PDF for m='+str(m) )
    pl.plot(results[0].keys(), results[0].values(), 'r-')
    pl.plot(results[1].keys(), results[1].values(), 'b-')
    pl.plot(results[2].keys(), results[2].values(), 'g-')
    pl.plot(results[3].keys(), results[3].values(), 'm-')
    pl.legend(('m=1', 'm=2', 'm=4', 'm=8'),loc='upper right')
    pl.savefig('degree_dist_PRA'+ext)
    pl.clf()

    f.close()
    f1.close()


"""-----PHASE 3, RANDOM WALK DEGREE DISTRIBUTION"""


#compare with previous deg dists
def phase3():
    base=1.1
    nets=[]
    results=[]
    m=3
    L=0
    it=10000
    av=10
    for i in range(6):
        res={}
        nets.append(BAmodel(m,True))
        for j in range(av):
            interes={}
            nets[i]=BAmodel(m,True)
            nets[i].walk_evolve(it,L)
            values=sorted( nw.degree(nets[i].G).values())
            interes=degree_freqdist(values  )
            for degree, freq in interes.items():
                if (degree in res):
                    res[degree]+=freq
                else:
                    res[degree]=freq
            print '%dth of 6 Graphs %d %% done' % (i+1,(j+1))
        for degree in res.keys():
            res[degree]/=av
        results.append(col.OrderedDict(sorted( f_to_p(res, it+m).items() )))
        L+=1
        if (i==4):
            L=6
        if (i==5):
            L=10

    pl.loglog()
    pl.ylabel('p(k)')
    pl.xlabel('k')
    pl.title('Degree PDF for m='+str(m) )
    pl.plot(results[0].keys(), results[0].values(), 'r-')
    pl.plot(results[1].keys(), results[1].values(), 'b-')
    pl.plot(results[2].keys(), results[2].values(), 'g-')
    pl.plot(results[3].keys(), results[3].values(), 'm-')
    pl.plot(results[4].keys(), results[4].values(), 'k-')
    pl.plot(results[5].keys(), results[5].values(), 'c-')
    pl.legend(('L=0', 'L=1', 'L=2', 'L=3','L=6','L=10'),loc='upper right')
    pl.savefig('degree_dist_RW'+ext)
    pl.clf()






"""----------------------------------MAIN PROGRAM-----------------------------------"""

#phase13()
phase14()
#phase2()
#phase3()




