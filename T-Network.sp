*** T network***
.temp27

***Supply voltages***
Vin 1 0 pulse(0 1.8 0s 100p 100p 5n 10n)

***network***
R 1 2 60
L 2 3 5n
C 2 0 120f
Cl 3 0 10f

**operation part**
.option post
.tran 1p 100n

**measure**
.measure tran tr trig v(1) val=0.9 rise=1 targ V(3) val=0.9 rise=1 
.measure tran tf trig v(1) val=0.9 fall=1 targ V(3) val=0.9 fall=1 
*.measure tran td  val='(tr+tf)/2'  

.measure tran Pavg AVG par('V(1)*I(Vin)') from 0 to 100n

**plot**
.print V(3) I(L) I(C)

**end**
.end

