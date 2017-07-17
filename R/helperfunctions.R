##Renske van Raaphorst
##7/14/2017

##Helperfunctions: functions necessary to make other functions working properly.

##Set the pixel to um conversion
spotrSetpixel2um <- function(){
  pixel2um <- readline(caption="Give the conversion factor from pixel to um./n")
  return(pixel2um)
}

##error message + solution when pixel2um is not set
spotrPixelerror <- function(){
  errormessage <- readline(caption="The conversion factor from pixel to um is not indicated. Please use 'spotrSetpixel2um' by pressing 'a'.
          if you don't want to convert from pixel to um, press 'b'./n")
  if(errormessage=="a"|errormessage=="A"){
    conv == pixel2um()
  }
  if(errormessage=="b"|errormessage=="B"){
    conv == 1
  }
  else(print("Did not receive 'a' nor 'b'. If you want to convert from pixel to um, use the function 'spotrSetpixel2um()' manually./n"))
  return(conv)
}

