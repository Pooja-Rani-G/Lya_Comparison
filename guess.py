import math

mem_per_node = 187
L = input('Enter level: ')

L = int(L)

ngridtot = 18*8**(L-2)
nparttot = 2*8**L

Mem_tot = 1.4*ngridtot/1e6 + 0.7 * nparttot/1e7

print('\n For a hydro-cosmo RAMSES simulation on DTP cluster ...\n')
print('ngridtot =',ngridtot)
print('nparttot =',nparttot)

print('Total memory required (GB) =',Mem_tot)
print('Number of nodes required = ',math.ceil(Mem_tot/mem_per_node))
