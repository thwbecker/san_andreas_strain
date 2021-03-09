#
# average values given as x y sigma_y, already sorted by x, with dx spacing
#
BEGIN{
    if(dx == 0)
	dx = 5;
    
}
{
    n++;
    x[n] = $1;
    y[n] = $2;
    s[n] = $3;
}
END{
    xmin = int((x[1]-dx/2)/dx) * dx;
    xmax = int((x[n]+dx/2)/dx) * dx + 1e-7;

    m = (xmax-xmin)/dx+1;
    for(i=1;i<=m;i++){
	a[i] = 0;
	as[i] = 0;
	an[i] = 0;
    }
    for(i=1;i <= n;i++){
	bin = 1 + int((x[i]-xmin)/dx);
	if(bin > m)print("error") > "/dev/stderr";
	a[bin] += y[i];
	as[bin] += s[i]**2;
	an[bin]++;
    }

    for(i=1;i <= m;i++)
	if(an[i]){
#	    print(xmin + (i-0.5)*dx, a[i]/an[i],sqrt(as[i])/an[i]);
	    print(xmin + (i-0.5)*dx, a[i]/an[i],sqrt(as[i])/sqrt(an[i]));
	}
}