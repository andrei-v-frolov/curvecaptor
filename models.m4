divert(-1) $Id: models.m4,v 1.8 2005/05/11 04:49:16 afrolov Exp $

# Curve Captor - vacuum tube curve capture and model builder tool
# Spice 3F4 model definitions - implemented as M4 macros
# 
# Copyright (C) 2001-2005 Andrei Frolov <frolov@cita.utoronto.ca>
# Distributed under the terms of GNU Public License.

# Spice syntax
define(`macro', `Spice 3F4')

# Vacuum diode models
define(`diode',		`Bp  P K  I=($1m)*uramp(V(P,K))^1.5')
define(`diode_cp',	`Bp  P K  I=($1m)*uramp(V(P,K)+($2))^1.5')

# Vacuum triode models
define(`triode',	`Bp  P K  I=($1m)*uramp(($2)*V(G,K)+V(P,K))^1.5')
define(`triode_cp',	`Bp  P K  I=($1m)*uramp(($2)*V(G,K)+V(P,K)+($3))^1.5')
define(`rydel4',	`Bp  P K  I=(($1m)+($2m)*V(G,K))*uramp(($3)*V(G,K)+V(P,K)+($4))^1.5')
define(`rydel5',	`Bp  P K  I=(($1m)+($2m)*V(G,K))*uramp(($3)*V(G,K)+V(P,K)+($4))^1.5 * V(P,K)/(V(P,K)+($5))')
define(`koren4',	`Bp  P K  I=($1m)*uramp(V(P,K)*ln(1.0+exp(($2)+($2)*($3)*V(G,K)/V(P,K)))/($2))^($4)')
define(`koren5',	`Bp  P K  I=($1m)*uramp(V(P,K)*ln(1.0+exp(($2)+($2)*($3)*V(G,K)/sqrt(($4k)+V(P,K)^2)))/($2))^($5)')
define(`koren6',	`Bp  P K  I=($1m)*uramp(V(P,K)*ln(1.0+($2)+exp(($3)+($3)*(($4)+($5m)*V(G,K))*V(G,K)/V(P,K)))/($3))^($6)')
define(`koren8',	`Bp  P K  I=($1m)*uramp(V(P,K)*ln(1.0+($2)+exp(($3)+($3)*(($4)+($5m)*V(G,K))*V(G,K)/sqrt(($6)^2+(V(P,K)-($7))^2)))/($3))^($8)')

# Inter-element capacitance
define(`cap', `    Cgp G P  $1
    Ci  G K  $2
    Co  P K  $3')

divert(0)dnl
