#/usr/bin/perl/

$R=0.3548387097;
for my $i (0..63)
{
	$SNR = $i/4;
	$sigma= 1 / (sqrt(2*$R*10**($SNR/10)));
	print "    $sigma, // ($i) SNR = $SNR, RATE = $R\n";
}
