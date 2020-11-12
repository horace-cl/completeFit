# completeFit
Implementation of the complete fit PDF(theta, mass; params)

PDF = Y_s S(\theta)*efficiency(\theta) S(\mass) + Y_b B(\theta) B(\mass)

<br>
Where 
<br>
    Y_s and Y_b            -   are signal and background yields respectivelly
<br>
    B(\theta) and B(\mass) -   are the shapes of the bakground on the mass and angular variables
<br>
    S(\mass)               -   are the shapes of the signal on the mass and angular variables
<br>
    S(\theta)              -   is given by the theory
<br>
    efficiency(\theta)     -   is the efficiency on the angular variable
    
 <br><br>

The main script is completeSymfit.py it implements the complete fit using pdf's defined on the othe scripts.

FitSym.py is WIP, it is intended to poduce several fits on one single script, using the output as input for next iterations. 
Argparse is being used here
