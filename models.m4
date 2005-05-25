divert(-1) $Id: models.m4,v 1.9 2005/05/25 00:56:29 afrolov Exp $

# Curve Captor - vacuum tube curve capture and model builder tool
# Spice 3F4 model definitions - implemented as M4 macros
# 
# Copyright (C) 2001-2005 Andrei Frolov <frolov@cita.utoronto.ca>
# Distributed under the terms of GNU Public License.

# Spice syntax
define(`macro', `Spice 3F4')

define(`Vp', `V(P,K)')
define(`Vg', `V(G,K)')
define(`Ip', `Bp  P K  I')

define(`pow2', `($1)^2')
define(`ppow', `uramp($1)^$2')

# Vacuum diode models
define(`diode',		`Ip=($1m)*ppow(Vp,1.5)')
define(`diode_cp',	`Ip=($1m)*ppow(Vp+($2),1.5)')

# Vacuum triode models
define(`triode',	`Ip=($1m)*ppow(($2)*Vg+Vp,1.5)')
define(`triode_cp',	`Ip=($1m)*ppow(($2)*Vg+Vp+($3),1.5)')
define(`rydel4',	`Ip=(($1m)+($2m)*Vg)*ppow(($3)*Vg+Vp+($4),1.5)')
define(`rydel5',	`Ip=(($1m)+($2m)*Vg)*ppow(($3)*Vg+Vp+($4),1.5) * Vp/(Vp+($5))')
define(`koren4',	`Ip=($1m)*ppow(Vp*ln(1.0+exp(($2)+($2)*($3)*Vg/Vp))/($2),($4))')
define(`koren5',	`Ip=($1m)*ppow(Vp*ln(1.0+exp(($2)+($2)*($3)*Vg/sqrt(($4k)+pow2(Vp))))/($2),($5))')
define(`koren6',	`Ip=($1m)*ppow(Vp*ln(1.0+($2)+exp(($3)+($3)*(($4)+($5m)*Vg)*Vg/Vp))/($3),($6))')
define(`koren8',	`Ip=($1m)*ppow(Vp*ln(1.0+($2)+exp(($3)+($3)*(($4)+($5m)*Vg)*Vg/sqrt(pow2($6)+pow2(Vp-($7)))))/($3),($8))')

# Inter-element capacitance
define(`cap', `dnl
    Cgp G P  $1
    Ci  G K  $2
    Co  P K  $3')

divert(0)dnl
