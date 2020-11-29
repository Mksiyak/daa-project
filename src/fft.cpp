#include <bits/stdc++.h>
using namespace std;
#define CP complex<double>
#define VCP vector<complex<double> >
vector <CP> universal; 
const double PI=3.1415926535897932; 
const int MAXN = (1<<15); 
const double eps=1e-9; 
double rounder (double x)
{
	if (x-floor(x)<eps)
	{
		return floor(x); 
	}
	return x;
}
VCP fft (VCP coefficients, vector <CP> eval)
{
	if (coefficients.size()==1)
	{
		VCP basecase; 
		basecase.push_back(coefficients[0]);
		return basecase; 
	}
	VCP coeff1, coeff2; 
	vector <CP> eval1, eval2;
	for (int g=0; g<coefficients.size(); g+=2) coeff1.push_back(coefficients[g]); 
	for (int g=1; g<coefficients.size(); g+=2) coeff2.push_back(coefficients[g]); 
	for (int g=0; g<eval.size()/2; g++) eval1.push_back(eval[g]*eval[g]);
	for (int g=eval.size()/2; g<eval.size(); g++) eval2.push_back(eval[g]*eval[g]);
	VCP first=fft(coeff1, eval1), second=fft(coeff2, eval2);
	VCP toreturn; 
	for (int g=0; g<coefficients.size(); g++)
	{
		if (g<coefficients.size()/2)
		{
			toreturn.push_back(first[g]+eval[g]*second[g]);
		}
		else
		{
			toreturn.push_back(first[g-coefficients.size()/2]+eval[g]*second[g-coefficients.size()/2]); 
		}
	}
	return toreturn; 
} 
vector <int> interpolate (VCP a, VCP b)
{
	VCP first=fft(a,universal);
	VCP second=fft(b,universal); 
	VCP third;
	for (int g=0; g<first.size(); g++) third.push_back(first[g]*second[g]);
	for (int g=0; g<universal.size(); g++) universal[g]=complex<double>(1,0)/universal[g]; 
	//cout << second[0] << '\n';
	VCP fincoff=fft(third, universal); 
	for (int g=0; g<universal.size(); g++) universal[g]=complex<double>(1,0)/universal[g]; 
	//for (int g=0; g<fincoff.size(); g++) cout << fincoff[g] << '\n'; cout<<'\n';
	vector <int> answer;
	for (int g=0; g<fincoff.size(); g++)
	{
		answer.push_back(real(fincoff[g]+0.5)/MAXN); 
	}	
	return answer; 
}
void fill (VCP & a)
{
	int m=a.size(); 
	for (int g=0; g<MAXN-m; g++) a.push_back(complex<double>(0,0));
}
int main()
{
	ios_base::sync_with_stdio(0); 
    freopen("input.txt","r",stdin);
	int T = 1;
	for (int g=0; g<T; g++)
	{
		universal.clear(); 
		for (int g=0; g<MAXN; g++)
		{
		universal.push_back(complex <double> (rounder(cos((double)g*(2*PI)/(double) MAXN)), rounder(sin((double)g*(2*PI)/(double) MAXN))));
		}
		VCP a, b; 
		int N; cin >> N; 
		for (int g=0; g<N+1; g++)
		{
			int t; cin >> t; a.push_back(complex<double> (t, 0)); 
		}
		for (int g=0; g<N+1; g++)
		{
			int t; cin >> t; b.push_back(complex<double> (t, 0)); 
		}
		reverse(a.begin(), a.end()); 
		reverse(b.begin(), b.end()); 
		fill (a); 
		fill(b);
		vector <int> ans=interpolate(a, b); 
		reverse(ans.begin(), ans.begin()+2*N+1);
		for (int g=0; g<2*N+1; g++) cout << ans[g] << ' '; 
		cout << '\n';
	}
	return 0; 
}
