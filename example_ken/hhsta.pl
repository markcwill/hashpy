: # use perl
eval 'exec $ANTELOPE/bin/perl -S $0 "$@"'
if 0;

use lib "$ENV{ANTELOPE}/data/perl" ;
use Datascope;
#use Math::Trig;
#use List::Util 'max';
#use List::Util 'min';
#use POSIX qw(log10);

$numargs = @ARGV;

if ($numargs < 1) {
   print ("\n");
        print "Usage: hhsta.pl [HASH Station File]\n";
        print "       create station file for HASH from dbmaster\n";
   die   ("\n");
  }

   open(STAFILE,">$ARGV[0]");
   $filename = $ARGV[0];

   $DB = "/home/ymp24/ANTELOPE/RENO/dbmaster/reno";

      if(-e $DB) {

# Join tables in dbprocess

      @db  = dbopen("$DB", "r");

      @stachans = dbprocess ( @db, "dbopen site",
                 "dbjoin sitechan",
                 "dbjoin affiliation");

       $numchans = dbquery( @stachans , "dbRECORD_COUNT" ); 
   
       print "Number of Sta Chans:  $numchans\n";

# Event summary information. 

     $stas = 0;

      while ($stas < $numchans) {
        $stachans[3] = $stas; 
        ($sta,$net,$chan,$lat,$lon,$elev,$staname,$ondate,$offdate) = dbgetv(@stachans,"site.sta","affiliation.net","sitechan.chan","site.lat","site.lon","site.elev","site.staname","site.ondate","site.offdate");
         &format_text;
#        print         "$sta $chan $staname $lat $lon$elev $ondate $offdate\n";
#        print STAFILE "$sta $chan $staname $lat $lon$elev $ondate $offdate\n";
         print         "$numchans-$stas $outputformat $onday $offday $net\n";
         print STAFILE "$outputformat $onday $offday $net\n";
         $stas++;
       }
        dbclose (@db);
      }  
        close(STAFILE);
        system("cat $filename | sort | uniq > stationfile.hash");

#    }

#--------------------------------------------------------------

   sub format_text {

         $elev = $elev * 1000;
#        $chan = sprintf("%-3.3s",$chan);
#        $staname = sprintf("%-32.32s",$staname);
#        $lat = sprintf("%8.5f",$lat);
#        $lon = sprintf("%10.5f",$lon);
#        $elev = sprintf("%4.4d",$elev);
#
         $onday  = `epoch +"%Y/%m/%d" $ondate`;
         if ($offdate < 0) { 
          $offdate = '2050001';
         }
          $offday  = `epoch +"%Y/%m/%d" $offdate`;
          chop($onday);
          chop($offday);
         $outputformat = sprintf("%-4.4s %-3.3s %-32.32s %8.5f %10.5f  %4.4d",$sta,$chan,$staname,$lat,$lon,$elev);
       }


#BBR  VHZ Big Bear Solar Observatory       34.26230 -116.92075  2037 2000/06/14 3000/01/01 CI
#BBR  VHZ Big Bear Solar Observatory       34.26230 -116.92075  2069 2000/06/14 3000/01/01 CI      
#BBR  VHZ SCSN STATION  134                34.24170 -116.90914  2016 1900/01/01 3000/01/01 CI
#BBRV VHZ SCSN STATION  134                34.24170 -116.90914  2016 1900/01/01 3000/01/01 CI
#BC2  EHZ BIG CHUCKAWALLA MTNS             33.65702 -115.46202   995 1977/01/01 1991/05/08 CI
#BC2  EHZ BIG CHUCKAWALLA MTNS             33.65702 -115.46202   995 1977/01/01 1991/05/08 CI      
#BCC  LLE Bear Creek Country Club          33.57508 -117.26119   363 2000/07/11 3000/01/01 CI
#BCC  LLE Bear Creek Country Club          33.57508 -117.26119   391 2000/07/11 3000/01/01 CI      
#BCC  LLN Bear Creek Country Club          33.57508 -117.26119   363 2000/07/11 3000/01/01 CI
#BCC  LLN Bear Creek Country Club          33.57508 -117.26119   391 2000/07/11 3000/01/01 CI      
