function sym_vc=add_vc(symbol,pm)
sym_vc(2:27)=[symbol(25:30) pm symbol(31:43) -pm symbol(44:48)];
sym_vc(39:64)=[symbol(1:5) pm symbol(6:18) pm symbol(19:24)];