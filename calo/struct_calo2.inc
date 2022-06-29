! STRUCT_CALO.INC

! Define some record structures

! ... where do we loose events?
	integer*4	calo2STOP_trials,calo2STOP_slit_hor,calo2STOP_slit_vert,calo2STOP_successes

! ... a general common block
	common /spectrometer_calo2/  calo2STOP_trials,calo2STOP_slit_hor
     > ,calo2STOP_slit_vert,calo2STOP_successes

