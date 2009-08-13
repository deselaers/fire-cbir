#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use File::Basename;
use File::Copy;
use LWP::Simple qw(!head);
use Cwd;
use File::Spec;

# absolute path to current script (auto)
my $script_dir=cwd;

# !!!!! EDIT HERE !!!!!
# absolute path to the WEB-server (EDIT manually)
my $web_server_dir='/home/demo/apache/htdocs';
# absolute path to the IMAGE-directory (EDIT manually)
my $pic_dir='/home/demo/apache/htdocs/img';
# absolute path to the HTML-directory (EDIT manually), which contains:
# select_rectangle.html
# js/select_rectangle.js
# js/wz_jsgraphics.js
my $html_dir='/home/demo/apache/htdocs/select_rectangle';
# !!!!! END EDIT !!!!!

# relative path from current script to HTML-directory (auto)
my $html_from_script_dir='../'.substr($html_dir, length($web_server_dir)+1);
# realtive path from current script to IMAGE-directory (auto)
my $pic_from_script_dir='../'.substr($pic_dir, length($web_server_dir)+1);

# forming current time format YYYY_MM_DD__HH_MM_SS
my @t=localtime;
my $time=(1900+$t[5]).'_'.sprintf("%02d", 1+$t[4]).'_'.sprintf("%02d", $t[3]).'__'.sprintf("%02d", $t[2]).'_'.sprintf("%02d", $t[1]).'_'.sprintf("%02d", $t[0]);
# function to insert suffix (is the first argument) before the file extension for the file-array (second argument)
sub add_to_name() {
 my $add=shift;
 my @names=@_;
 for(@names) { s!(\.[^\.]+)$!__${add}$1!; }
 if($#names>0) { return @names; } else { return $names[0]; }
}

# error detection, default empty string means no errors
my $error='';

# printing the error log message (first argument - message, second - reason)
# time of execution is added to the begin of the message
sub error_log() {
 my ($message, $reason)=@_;
 # error log file. If error occurs, the exact reason will be written
 open ERROR_LOG, ">>${html_dir}/error_log.txt";
 print ERROR_LOG ("$time:\t$message\n".(' ' x 21)."\tReason: $reason\n");
 close(ERROR_LOG);
 $error="ERROR:\n$message\n\nREASON:\n$reason\n";
 &ReturnHTML;
}

# maximum size limit for the file to upload in Kbytes
my $max_file_size=5 * 1024; # 5 Mbytes
# additional internal POST protection
$CGI::POST_MAX=($max_file_size + 40) * 1024;

# initialization of CGI
my $query = new CGI;
my ($param, $to_do, $upload_type, $port, $image_name, $image_name_perm, $cksum, $image_name_selected)=('', '', '', '', '', '', '');

# without arguments, script produces the start page
unless($query->param()) {
 if(cgi_error) {
  &error_log("CGI external error.\n(Perhaps maximum file size (${max_file_size} Kb) is exceeded)", cgi_error);
 }
 else {
  &ReturnHTML;
 }
}

# PARAMETERS
# to_do = 'UpdateImage' or 'ShowSelectedPart_and_SaveCoordinates'
$param=$query->param('to_do');
if(defined $param) { $to_do=$param; }
# verification of the 'to_to' parameter
if( $to_do!~m!^(UpdateImage|ShowSelectedPart_and_SaveCoordinates)$!i) {
 &error_log("Cannot proceed. Internal server (or browser) error.", "CGI-argument 'to_do' is '${to_do}'.\nBut it has to be 'UpdateImage' or 'ShowSelectedPart_and_SaveCoordinates'");
}
# upload_type = 'www', 'local' or 'server'
$param=$query->param('upload_type');
if(defined $param) { $upload_type=$param; }
# verification of the 'upload_type' parameter
if( $upload_type!~m!^(www|local|server)$!i) {
 &error_log("Cannot proceed. Internal server (or browser) error.", "CGI-argument 'upload_type' is '${upload_type}'.\nBut it has to be 'www', 'local' or 'server'.\nTry to check the radio-button in '1) Load image...'-section.");
}
# port = port for server
$param=$query->param('port');
if(defined $param) { $port=$param; }
my $port_default='80';
# check sum verification
if($port!~m!^\d+$!) {
 $port=$port_default;
}
# image source file name
$param=$query->param('image_name');
if(defined $param) { $image_name=$param; }
# security: as input image name, one may use only letters, numbers, underscores, spaces, periods, exclamation, points, question marks, and hyphens.
my $image_name_OK=1;
if($upload_type eq 'local') {
 if($image_name=~m!^[\w\\\/ .:\!?-]+$!) {
  if( ($image_name=~m!\/!) && ($image_name=~m!\\!) ) {
   $image_name_OK=0;
  }
  else {
  $image_name=~s!\\!\\\\!g;
  }
 }
 else {
  $image_name_OK=0;
 }

 if($image_name_OK != 1) {
  my $image_name_err=$image_name;
  $image_name_err=~s![\w .:\!?-]+!!g;
  if($image_name_err=~m!\/!) {
   $image_name_err=~s!\/!!g;
  }
  else {
   $image_name_err=~s!\\!!g;
  }
  &error_log("File \"$image_name\" was not copied to server.", ("File name includes (an) invalid character(s):\n".$image_name_err."\nYou may only enter letters, numbers, underscores, spaces, periods, exclamation, points, question marks, and hyphens.") );
 }
}
else {
$image_name_OK=1;
}

# main image file name on server ('perm' = permanent)
$param=$query->param('image_name_perm');
if(defined $param) { $image_name_perm=$param; }
if( ($to_do eq 'ShowSelectedPart_and_SaveCoordinates') && (!-e $image_name_perm) ) {
 &error_log("Cannot crop the selection. Internal server error.", "Image does not exist on server");
}
# check sum of the main image file (UNIX command: cksum <file_name>)
$param=$query->param('cksum');
if(defined $param) { $cksum=$param; }
my $cksum_default='not_defined';
# check sum verification
if($cksum!~m!^\d+_\d+$!) {
 $cksum=$cksum_default;
}
# directory, where the pictures have to be saved
# my $image_dir='pic/';
my $image_dir=$pic_dir.'/';
if(! -e $image_dir) { mkdir($image_dir); }

my $basename=basename($image_name);
# check for MS Windows file name
if($image_name=~m!\\([^\\]+)$!) {
  $basename=$2;
}
# temporary main image file name on server
my $image_name_tmp=$image_dir.$basename;

# crop image file name on server
$image_name_selected=&add_to_name($cksum.'__'.$time, $image_name_tmp);

# MAIN processing
if($to_do eq 'UpdateImage') {
  # uploading main image file to server
  if($upload_type eq 'www') {
   my @head=head($image_name); # head format: ($content_type, $document_length, $modified_time, $expires, $server)
   if($head[1]>$max_file_size*1024) {
    &error_log("File \"$image_name\" (WWW-link) was not copied to server.", "Maximum file size (${max_file_size} Kb) is exceeded");
   }
   my $http_response_code=getstore($image_name, $image_name_tmp);
   if(! -e $image_name_tmp) {
    &error_log("File \"$image_name\" (WWW-link) was not copied to server.", "HTTP response code: \"$http_response_code\"");
   }
  }
  elsif($upload_type eq 'local') {
   my $file = $query->param('file_name_local');
   my $size=0;
   open DAT,'>'.$image_name_tmp or &error_log("Internal server error: during processing temporary file.", $!);
   binmode $file;
   binmode DAT;
   my $data;
   while(read $file,$data,1024) {
     $size++;
     if ($size>$max_file_size) {
      close DAT;
      unlink $image_name_tmp;
      &error_log("Local file \"$image_name\" was not copied to server.", "Maximum file size (${max_file_size} Kb) is exceeded");
      last;
     }
     print DAT $data;
   }
   close DAT;
  }
  else {
   copy($image_name, $image_name_tmp);
  }
  # existence of the image
  if(! -e $image_name_tmp) {
   &error_log("File \"$image_name\" was not copied (or found) to server.", "Unknown internal server error");
  }
  # detection of the file type
  chomp(my $file_type=`file -b "$image_name_tmp"`);
  if($file_type=~/(image|bitmap)\sdata/i) {
   # computing the CheckSum of the copied file
   $cksum=`cksum "$image_name_tmp"`;
   if($cksum=~m!^(\d+)\s(\d+).*$!) {
    $cksum=$1.'_'.$2;
   }
   else {
    $cksum=$cksum_default;
   }
   # forming the permanent main image file name on server (including CheckSum and Time)
   $image_name_perm=&add_to_name($cksum.'__'.$time, $image_name_tmp);
   rename($image_name_tmp, $image_name_perm);
  }
  elsif(-s $image_name_tmp == 0) {
   unlink $image_name_tmp;
   &error_log("File \"$image_name\" was not copied to server.", "File is empty or does not exist");
  }
  else {
   unlink $image_name_tmp;
   &error_log("File \"$image_name\" was not copied to server.\nIts type is \"${file_type}\".", "Only images can be copied");
  }
}
else {
 # to_do = 'ShowSelectedPart_and_SaveCoordinates'
 my ($coord, %c)=('');
 # forming the coordinates TXT-file, all agruments of this CGI-script will be saved, including coordinates
 # arguments are also saved in hash-array %c
 open COORD, '>'.(&add_to_name($cksum.'__'.$time, $image_name_tmp).'.txt');
 for $coord($query->param) {
  $c{$coord}=$query->param($coord);
  print COORD ("$coord=", $c{$coord}, "\n");
 }
 close(COORD);

 # converting the main image file name on server (croping the selected rectangle)
 unlink $image_name_selected;
 system("convert \"$image_name_perm\" -crop ".($c{'x_end'}-$c{'x_start'}).'x'.($c{'y_end'}-$c{'y_start'}).'+'.$c{'x_start'}.'+'.$c{'y_start'}." \"$image_name_selected\"");
}

&ReturnHTML;

sub ReturnHTML {
# HTML document, which will be returned to browser
print header;
my $html_source='select_rectangle.html';
if($upload_type eq 'server') {
 $html_source='select_rectangle_fire.html';
}
if(! open MAIN_HTML, "<${html_dir}/${html_source}") {
 print h2("ERROR:<br>HTML-file<br>${html_from_script_dir}/${html_source} <br> was not found on server");
}
my ($image_name_perm_browser, $image_name_selected_browser)=($image_name_perm, $image_name_selected);
for ($image_name_perm_browser, $image_name_selected_browser) {
 s/${pic_dir}/${pic_from_script_dir}/;
}

my $line;
for $line(<MAIN_HTML>) {
 $line=~s!(type="text\/javascript" src=\")(js\/)!$1${html_from_script_dir}\/$2!; # $line=~s!(form action=\")(save_coord.cgi\")!$1${script_from_html_dir}\/$2!;
 if( ($line=~s!(\<img id="main_image")!$1 src="${image_name_perm_browser}"!) &&
     ($to_do eq 'ShowSelectedPart_and_SaveCoordinates') )
  {
    $line.="<br>Selected part:<br><img src=\"${image_name_selected_browser}\">";    $line.="<form method=\"post\" name=\"queryForm\" action=\"fire-parts.cgi\"> <input type=\"hidden\" name=\"serverfile\" value=\"$image_name_selected\"> <input type=\"hidden\" name=\"port\" value=\"$port\"> <input type=\"submit\" value=\"search similar images\"> </form>";

  }
 # recovering the cheched image source radio button on webpage ('www' or 'local')
 $line=~s!(name="upload_type" value="${upload_type}")!$1 checked="checked"!;
 # main initialization of the JavaScript
 if(  ($line=~s!(InitMain\(')(',')(',')(',')(',')(',')('\);)!$1${upload_type}$2${image_name}$3${image_name_perm}$4${cksum}$5${image_name_selected}$6${port}$7!) &&
 ($error ne '') ) {
  $error=~s!([^\w\d\n])!\\$1!g;
#  $error=~s!\\!\\\\!g;
#  $error=~s!"!\\"!g;
  $error=~s!\n!\\n!g;
  $line.="alert(\"${error}\");\n";
 }

 # recovering the value of the image source field.
 # NOTE: Browsers are typically ignoring the value for the FileSelect object, due to security reasons
 $line=~s!(\<input id="file_name_${upload_type}")!$1 value="$image_name"!;
 print $line;
}
close(MAIN_HTML);
exit;
}
