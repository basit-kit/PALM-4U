	  �  B   k820309              15.0        ֳ�X                                                                                                           
       smog_Rates.f90 SMOG_RATES                                                    
                                                          
                                                                                                                                                             
                                                       
                                                       
                  @                                    
                  @                                                  
      p          p            p                          %         @                                	                   
       #ARR%EXP 
   #ARR%DBLE    #A0    #B0    #C0                                              
     EXP                                                DBLE            @                                    	                  @                                    	                  @                                    	       %         @                                                   
       #ARR2%EXP    #ARR2%DBLE    #A0    #B0                                                   EXP                                                DBLE            @                                    	                  @                                    	       %         @                                                   
       #EP2%EXP    #EP2%DBLE    #A0    #C0    #A2    #C2    #A3    #C3                                                   EXP                                                DBLE            @                                    	                  @                                    	                  @                                    	                  @                                    	                  @                                    	                  @                                    	       %         @                                                   
       #EP3%EXP    #EP3%DBLE    #A1     #C1 !   #A2 "   #C2 #                                                  EXP                                                DBLE            @                                     	                  @                               !     	                  @                               "     	                  @                               #     	       %         @                                $                   
       #FALL%LOG10 %   #FALL%EXP &   #FALL%DBLE '   #A0 (   #B0 )   #C0 *   #A1 +   #B1 ,   #C1 -   #CF .                                             %     LOG10                                           &     EXP                                           '     DBLE            @                               (     	                  @                               )     	                  @                               *     	                  @                               +     	                  @                               ,     	                  @                               -     	                  @                               .     	       %         H                               /                   
       #K_3RD%LOG10 0   #TEMP 1   #CAIR 2   #K0_300K 3   #N 4   #KINF_300K 5   #M 6   #FC 7                                             0     LOG10           
                                 1     
                
                                 2     
                
                                  3     	                
                                  4     	                
                                  5     	                
                                  6     	                
                                  7     	      %         H                               8                   
       #K_ARR%EXP 9   #K_298 :   #TDEP ;   #TEMP <                                             9     EXP           
                                  :     	                
                                  ;     	                
                                 <     
      #         @                                   =                    #UPDATE_SUN%COS >   #UPDATE_SUN%INT ?                                             >     COS                                           ?     INT #         @                                   @                     #         @                                   A                        �   "      fn#fn     �   @   J   SMOG_PARAMETERS      @   J   SMOG_GLOBAL "   B  p       DP+SMOG_PRECISION !   �  @       TEMP+SMOG_GLOBAL $   �  @       CFACTOR+SMOG_GLOBAL !   2  @       TIME+SMOG_GLOBAL     r  @       SUN+SMOG_GLOBAL #   �  �       RCONST+SMOG_GLOBAL    F  �       ARR    �  <      ARR%EXP      =      ARR%DBLE    B  @   a   ARR%A0    �  @   a   ARR%B0    �  @   a   ARR%C0      }       ARR2      <      ARR2%EXP    �  =      ARR2%DBLE    �  @   a   ARR2%A0    8  @   a   ARR2%B0    x  �       EP2      <      EP2%EXP    O  =      EP2%DBLE    �  @   a   EP2%A0    �  @   a   EP2%C0      @   a   EP2%A2    L  @   a   EP2%C2    �  @   a   EP2%A3    �  @   a   EP2%C3    	  �       EP3    �	  <      EP3%EXP    �	  =      EP3%DBLE    
  @   a   EP3%A1    P
  @   a   EP3%C1    �
  @   a   EP3%A2    �
  @   a   EP3%C2      �       FALL    �  >      FALL%LOG10      <      FALL%EXP    ?  =      FALL%DBLE    |  @   a   FALL%A0    �  @   a   FALL%B0    �  @   a   FALL%C0    <  @   a   FALL%A1    |  @   a   FALL%B1    �  @   a   FALL%C1    �  @   a   FALL%CF    <  �       K_3RD    �  >      K_3RD%LOG10    !  @   a   K_3RD%TEMP    a  @   a   K_3RD%CAIR    �  @   a   K_3RD%K0_300K    �  @   a   K_3RD%N     !  @   a   K_3RD%KINF_300K    a  @   a   K_3RD%M    �  @   a   K_3RD%FC    �  ~       K_ARR    _  <      K_ARR%EXP    �  @   a   K_ARR%K_298    �  @   a   K_ARR%TDEP      @   a   K_ARR%TEMP    [  p       UPDATE_SUN    �  <      UPDATE_SUN%COS      <      UPDATE_SUN%INT    C  H       UPDATE_RCONST    �  H       UPDATE_PHOTO 