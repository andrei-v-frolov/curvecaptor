divert(-1) $Id: models.m4,v 1.4 2001/09/20 04:26:52 frolov Exp $

# Curve Captor - vacuum tube curve capture and model builder tool
# Spice 3F4 model definitions - implemented as M4 macros
# 
# Copyright (C) 2001 Andrei Frolov <andrei@phys.ualberta.ca>
# Distributed under the terms of GNU Public License.

# Vacuum diode models
define(`diode',		`Bp  P K  I=($1m)*uramp(V(P,K))^1.5')
define(`diode_cp',	`Bp  P K  I=($1m)*uramp(V(P,K)+($2))^1.5')
define(`perugini',	`Bp  P K  I=(($1m)+($2m)*V(P,K))*uramp(V(P,K)+($3))^($4)')

# Vacuum triode models
define(`triode',	`Bp  P K  I=($1m)*uramp(($2)*V(G,K)+V(P,K))^1.5')
define(`triode_cp',	`Bp  P K  I=($1m)*uramp(($2)*V(G,K)+V(P,K)+($3))^1.5')
define(`rydel4',	`Bp  P K  I=(($1m)+($2m)*V(G,K))*uramp(($3)*V(G,K)+V(P,K)+($4))^1.5')
define(`rydel5',	`Bp  P K  I=(($1m)+($2m)*V(G,K))*uramp(($3)*V(G,K)+V(P,K)+($4))^1.5 * V(P,K)/(V(P,K)+($5))')
define(`rydelg',	`Bg  G K  I=($1m)*uramp(V(G,K))^1.5 * ((($2)+V(P,K))/(($3)+V(P,K)))^4')
define(`scott',		`Bp  P K  I=($1m)*uramp(ln(1.0+exp(($2)*(($3)*V(G,K)+V(P,K))))/($2))^($4)')
define(`koren4',	`Bp  P K  I=($1m)*uramp(V(P,K)*ln(1.0+exp(($2)+($2)*($3)*V(G,K)/V(P,K)))/($2))^($4)')
define(`koren5',	`Bp  P K  I=($1m)*uramp(V(P,K)*ln(1.0+exp(($2)+($2)*($4)*V(G,K)/sqrt(($3k)+V(P,K)^2)))/($2))^($5)')
define(`koren6',	`Bp  P K  I=($1m)*uramp(V(P,K)*ln(1.0+exp(($2)+($2)*(($4)+($5m)*V(G,K))*V(G,K)/V(P,K)))/($2)+($3)*V(G,K))^($6)')

# Inter-element capacitance
define(`cap', `    Cgp G P  $1
    Ci  G K  $2
    Co  P K  $3')

divert(0)dnl
