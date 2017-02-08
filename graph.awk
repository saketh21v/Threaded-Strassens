
BEGIN {	
	f1 = "recurGraph.tr"
	f2 = "strassenGraph.tr"
}

{

	num = $1;

	t1 = $2;

	t2 = $3;

	
	printf("%d	%lf\n", num, t1) > f1;
	printf("%d	%lf\n", num, t2) > f2;
	
}

END {
}
