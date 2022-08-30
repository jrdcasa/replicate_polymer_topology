import mbuild as mb

ethane_UA = mb.Compound()
ch3_1 = mb.Particle(name='_CH3', pos=[0, 0, 0])
ch3_2 = mb.Particle(name='_CH3', pos=[0.15, 0, 0])
ethane_UA.add([ch3_1, ch3_2])
ethane_UA.add_bond((ch3_1, ch3_2))

ethane_UA_box = mb.fill_box(ethane_UA, 100, box=[2, 2, 2])
ethane_UA_box.save('ethane-UA-box.gro', overwrite=True)
ethane_UA_box.save('ethane-UA-box.top', forcefield_name='trappe-ua', overwrite=True)