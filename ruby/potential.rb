
class Potential
  
  EV_TO_ERGS=1.60217646e-12
  
  def initialize(xmin, xmax, xstep, proc)
    @xdata = []
    xmin.step(xmax, xstep) {|x| 
      @xdata << x
      }
    @ydata = [1]
    @data = []
    @xdata.each{|x|
      @data << proc.call(x)*EV_TO_ERGS
    }
  end
  
  def xdata 
    @xdata
  end
  
  def write_binary(filename)
    format = "E*" # double precision little-endian
    File.open(filename, 'w') do |f|
      f.print([@xdata.size].pack(format))
      f.print([@ydata.size].pack(format))
      f.print(@xdata.pack(format))
      f.print(@ydata.pack(format))
      f.print(@data.pack(format))
    end
  end
  
end