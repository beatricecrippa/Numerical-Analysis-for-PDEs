# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 12:25:53 2020

@author: beacr
"""

from problems import *

dir = 'p3'

rblossFull = []
rlossFull = []
rl2Full = []

for layer in [1,2,3]:
		npde = ProblemBLSingularity_BD(64, layer, 50) # works very well
		with tf.Session() as sess:
			sess.run(npde.init)
			for i in range(20000):
				if( i>1000 ):
					break
				npde.train(sess, i)
			#		if i%100==0:
			#			npde.visualize(sess, False, i=i, savefig=dir)
				
		
		rblossFull.append(npde.rbloss)
		rlossFull.append(npde.rloss)
		rl2Full.append(npde.rl2)
				
plt.close('all')
plt.figure(1)
#plt.semilogy(npde.rbloss)
for i in range(len(rblossFull)):
        plt.semilogy(rblossFull[i],label="nLayer=%d" %layers[i])
plt.xlabel('Iteration')
plt.ylabel('$L_b$')
plt.legend()
plt.savefig(dir + str(dim) + '_' 'lb.png')
plt.savefig('./high.p1/lb.png')

#plt.close('all')
plt.figure(2)
#plt.semilogy(npde.rloss)
for i in range(len(rlossFull)):
        plt.semilogy(rlossFull[i],label="nLayer=%d" %layers[i])
plt.xlabel('Iteration')
plt.ylabel('$L_i$')
plt.legend()
plt.savefig(dir + str(dim) + '_' 'lb.png')
plt.savefig('./high.p1/li.png')

#plt.close('all')
plt.figure(3)
#plt.semilogy(npde.rl2)
for i in range(len(rl2Full)):
        plt.semilogy(rl2Full[i],label="nLayer=%d" %layers[i])
plt.xlabel('Iteration')
plt.ylabel('$||u-u_h||_2$')
plt.legend()
plt.savefig(dir + str(dim) + '_' 'lb.png')
plt.savefig('./high.p1/l2.png')