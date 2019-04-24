#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

my $mode = shift;
if(! $mode){
	my $src = basename($0);
	print STDERR "[Error] Incorrect command\n";
	print STDERR " Usage: \$ $src install or $src uninstall\n"; 

	exit;
}

my $install_dir = dirname($0);
$install_dir = abs_path($install_dir);

my $lastz_dir = "$install_dir/src/lastz-distrib-1.04.00/";
my $makeBlocks_dir = "$install_dir/src/makeBlocks/";
my $kent_dir = "$install_dir/src/kent/";
my $bedtools_dir = "$install_dir/src/bedtools2";

if($mode eq "install"){
	system "tar xvzf $install_dir/src.tar.gz";

	system "make -C $lastz_dir";
	system "make -C $makeBlocks_dir";
	system "make -C $bedtools_dir";
	
	system "cp path.conf.template path.conf";
	system "sed -i \"s|^lastz=\.*|lastz=$lastz_dir/src/lastz|\" $install_dir/path.conf";
	system "sed -i \"s|^makeBlocks=\.*|makeBlocks=$makeBlocks_dir|\" $install_dir/path.conf";
	system "sed -i \"s|^kent=\.*|kent=$kent_dir|\" $install_dir/path.conf";
	system "sed -i \"s|^bedtools=\.*|bedtools=$bedtools_dir/bin/bedtools|\" $install_dir/path.conf";
}elsif($mode eq "uninstall"){
	system "rm -f path.conf";

	system "make clean -C $lastz_dir";
	system "make clean -C $makeBlocks_dir";
	system "make clean -C $bedtools_dir";
	
	system "tar cvzf $install_dir/src.tar.gz";
	system "rm -rf src";
}else{
	my $src = basename($0);
	print STDERR "[Error] Incorrect command\n";
	print STDERR " Usage: \$ $src install or $src uninstall\n"; 

	exit;
}
