precipitation_data:
	* X consists of altitude, longitude, latitude (3)
	* Y consists of precipitation in Jan, Feb, ... Dec (12)
	* X -> Y authors renormalized X
	* Our result: (score 06.01.17) X->Y with epsilon = 0.008862
	
Chemnitz_data:
	* X consists of sin(phi_wind), cos(phi_wind), T (3)
	* Y consists of ozone, sulfur dioxid, dust, CO, NO_2, NO_x (7)
	* only use 1:1440 -- there are no more
	* X -> Y
	* Our result: (score 06.01.17) Result: X->Y with epsilon = 0.008506
	
Stock_returns:
	* X consists of the Stock returns SH,HSI,TWI,N225 (Asia 4 -- 1:4)
	* Y consists of the Stock retuns FTSE,DAX,CAC (Europe 3 -- 5:7) --> Stock_7
	* Y consists of the Stock retuns FTSE,DAX,CAC,DJ,NAS (Europe & USA 5 -- 5:9) --> Stock_9
	* X -> Y
	* Our result 7: Y->X with epsilon = 0.000558 (-r 0.001)
	* Our result 9: Y->X with epsilon = 0.000257 (-r 0.001)
	* ==> Normalizing data to fit standard normal:
	* 7: Result: X->Y with epsilon = 0.000050
	* 9: Result: X->Y with epsilon = 0.000532