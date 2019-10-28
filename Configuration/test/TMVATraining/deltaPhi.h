float deltaPhi(float phi1, float phi2) {
	float result=phi1-phi2;
	if (result>3.14) result-=6.28;
	if (result<=-3.14) result+=6.28;
	return result;
}
