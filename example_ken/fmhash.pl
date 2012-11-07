: # use perl
eval 'exec $ANTELOPE/bin/perl -S $0 "$@"'
if 0;

use lib "$ENV{ANTELOPE}/data/perl" ;
use Datascope;
use Math::Trig;
use List::Util 'max';
use List::Util 'min';
use POSIX qw(log10);

$numargs = @ARGV;

if ($numargs < 3) {
   print ("\n");
        print "Usage: fmhash.pl [Year] [Jday] [Orid]\n";
   die   ("\n");
  }

#Example '>fm.pl 2008 114 507701'     M 2.9 near Mogul 2008

   $year = $ARGV[0];
   $jday = $ARGV[1];
   $evorid = $ARGV[2];

   $jday = sprintf("%3.3d",$jday);
   $DB = "/data\/$year\/$jday/reno";

      if(-e $DB) {

# Join tables in dbprocess

      @db  = dbopen("$DB", "r");

      @dbday = dbprocess ( @db, "dbopen origin",
                 "dbsubset orid == $evorid",
                 "dbjoin assoc",
                 "dbjoin arrival",
                 "dbjoin site",
                 "dbjoin affiliation",
                 "dbsubset iphase =~ /.*[PS].*/");

# Count all phases
#     @dbsub = dbsubset( @dbday , "fm =~ /.*[cd].*/");

      $fms = dbquery( @dbday , "dbRECORD_COUNT" );   
      print " $DB $fms $evorid\n";

	  print "Phase Sta Chan Arid  Phase FM Def  Arr-Time     TT(sec)    Res     Dist  Az          SLat   SLon    Origin-Time          Elat   ELon       Depth  Mag  Database\n";
	  print "---------------------------------------------------------------------------------------------------------------------------------------------------------------\n";

      $ifms = 0; 

# Event summary information. 

      while ($ifms < $fms) {
        $dbday[3] = $ifms; 
        ($sta,$chan,$arid,$iphase,$phtime,$fm,$timeres,$ordtime,$jdate,$orid,$lat,$lon,$ml,$depth) = dbgetv(@dbday,"arrival.sta",
                "arrival.chan","arrival.arid","arrival.iphase","arrival.time","arrival.fm","assoc.timeres","origin.time","origin.jdate",
                "origin.orid","origin.lat","origin.lon","origin.ml","origin.depth");
        ($delta,$esaz,$timedef,$slat,$slon,$nass,$ndef) = dbgetv(@dbday,"assoc.delta","assoc.esaz","assoc.timedef","site.lat",
                "site.lon","origin.nass","origin.ndef");
         $traveltime = $phtime - $ordtime;
         $fm = substr($fm, 0, 1);
         $delta = $delta * 111.1;
         &format_text;
         print "$ifms $sta $chan $arid $iphase $fm $timedef $phtime $traveltime $timeres : $delta $esaz : $slat $slon : $ordtime $lat $lon $depth $ml $DB\n";
         $ifms++;
       }
        dbclose (@db);
      }  
#    }

#--------------------------------------------------------------

   sub format_text {

         ($latd,$latm) = unpack("A2 A5",$lat);
         ($lond,$lonm) = unpack("A4 A5",$lon);

         print "Source Location:  $latd $latm $lond $lonm\n";

         $latm = $latm * 60;
         $lonm = $lonm * 60;
         $latm = sprintf("%5.2f",$latm);
         $lonm = sprintf("%5.2f",$lonm);

         print "Source Loc Format:  $latd $latm $lond $lonm\n";

         $depth = sprintf("%05.2f",$depth);
         $timeres = sprintf("%7.4f",$timeres);
         $delta = sprintf("%6.2f",$delta);
         $sta = sprintf("%-4.4s",$sta);
         $esaz = sprintf("%3.3d",$esaz);
         $ifms = sprintf("%2.2d",$ifms);
         $lat = sprintf("%7.4f",$lat);
         $lon = sprintf("%9.4f",$lon);
         $slat = sprintf("%7.4f",$slat);
         $slon = sprintf("%9.4f",$slon);
         $traveltime = sprintf("%8.4f",$traveltime);
         $phtime   = `epoch +"%Y%m%e%H%M%S.%s" $phtime`;
         $ordtime  = `epoch +"%Y%m%e%H%M%S.%s" $ordtime`;
         chop($phtime);
         chop($ordtime);
         $phtime  = sprintf("%15.2f",$phtime);
         $ordtime = sprintf("%15.2f",$ordtime);
       }
