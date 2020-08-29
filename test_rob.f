      program circle
      real r, area
 
 
      write (*,*) 'Give Z r:'
      read  (*,*) r
      area = SQRT(0.77*0.77 +(0.93*0.93-0.77*0.77)/18*(r-18))
      write (*,*) 'RG = ', area
 
      stop
      end    
    