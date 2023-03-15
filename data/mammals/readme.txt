Data has been extracted from mammals data set. See description for the UICN_CONCERV_STATUS entries from www.iucnredlist.org for names.

X consists of environmental factors and Y is binary expressing the presence or absence of a species. The location information was removed from the data sets.

X:
elevation, precipitation, temperature, annual temperature range
(real valued)

Y (canis):
CANIS_AUREUS, CANIS_LUPUS
(binary)

Y (lepus):
LEPUS_EUROPAEUS, LEPUS_GRANATENSIS, LEPUS_TIMIDUS
(binary)

Y (martes)
MARTES_FOINA, MARTES_MARTES
(binary)

Y(canis_lepus_martes):
CANIS_AUREUS, CANIS_LUPUS, LEPUS_EUROPAEUS, LEPUS_GRANATENSIS, LEPUS_TIMIDUS, MARTES_FOINA, MARTES_MARTES
(binary)

Assumed ground truth for all of them:
X -> Y