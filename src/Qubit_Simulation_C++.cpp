#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
/*qubit中的複數的虛數部份以及實數部分都介於-1~1
 * 因此使用long long配合定點數可以有效的解決運算中的浮點誤差問題
 * 以2^62為縮放基準
 */

class Qubit_Simulation{

	public:

		static constexpr const long long int Fixed_Point=(1LL<<62);
		static constexpr const int Fixed_shift=62;
		static constexpr const double Root_half=(1.0/sqrt(2.0));
		static constexpr const double Pi=3.141592658979;
		static constexpr const unsigned long long int Fixed_Pi=14488038916154245120ULL;

		static long long int FixedPoint(const double input){

			return static_cast<long long>(input*Fixed_Point);

		}

		static long long int NewtonSqrt(__int128_t input){

			if(input==0)
				return 0;

			__int128_t Result=input;
			__int128_t x=(Result+1)/2;

			while(x<Result){

				Result=x;
				x=(Result+input/Result)/2;

			}

			return static_cast<long long int>(Result);

		}

        static double FixedPointToDouble(const long long int input){

			double result=1.0*input/Fixed_Point;

			return result;

		}

	private:

	class complex{

		public:

		long long int Real,Imaginary;

			complex():Real(0),Imaginary(0){};
			explicit complex(double A):Real(FixedPoint(A)),Imaginary(0){};
			explicit complex(double A,double B):Real(FixedPoint(A)),Imaginary(FixedPoint(B)){};
			complex(long long int A):Real(A),Imaginary(0){};
			complex(long long int A,long long int B):Real(A),Imaginary(B){};

			friend complex operator+(complex& c,const double d)noexcept{
				return complex(c.Real+FixedPoint(d),c.Imaginary);
			}
			friend complex operator+(const double d,complex& c)noexcept{
				return complex(c.Real+FixedPoint(d),c.Imaginary);
			}
			friend complex operator+(const complex& c,const long long int d)noexcept{
				return complex(c.Real+d,c.Imaginary);
			}
			friend complex operator+(const long long int d,const complex& c)noexcept{
				return complex(c.Real+d,c.Imaginary);
			}
			complex operator+(const complex& other)const& noexcept{
				return complex(Real+other.Real,Imaginary+other.Imaginary);
			}


			friend complex operator-(const complex& c,const double d)noexcept{
				return complex(c.Real-FixedPoint(d),c.Imaginary);
			}
			friend complex operator-(const double d,const complex& c)noexcept{
				return complex(FixedPoint(d)-c.Real,c.Imaginary);
			}
			friend complex operator-(const complex& c,const long long int d)noexcept{
				return complex(c.Real-d,c.Imaginary);
			}
			friend complex operator-(const long long int d,const complex& c)noexcept{
				return complex(d-c.Real,c.Imaginary);
			}
			complex operator-(const complex& other)const& noexcept{
				return complex(Real-other.Real,Imaginary-other.Imaginary);
			}


			friend complex operator*(const complex& c,const double d)noexcept{
				__int128_t Real_result=((__int128_t)c.Real*FixedPoint(d))>>Fixed_shift,
						   Imaginary_result=((__int128_t)c.Imaginary*FixedPoint(d))>>Fixed_shift;
				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend complex operator*(const double d,const complex& c)noexcept{
				__int128_t Real_result=((__int128_t)c.Real*FixedPoint(d))>>Fixed_shift,
						   Imaginary_result=((__int128_t)c.Imaginary*FixedPoint(d))>>Fixed_shift;
				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend complex operator*(const complex& c,const long long d)noexcept{
				__int128_t Real_result=((__int128_t)c.Real*d)>>Fixed_shift,
						   Imaginary_result=((__int128_t)c.Imaginary*d)>>Fixed_shift;
				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend complex operator*(const long long d,const complex& c)noexcept{
				__int128_t Real_result=((__int128_t)c.Real*d)>>Fixed_shift,
						   Imaginary_result=((__int128_t)c.Imaginary*d)>>Fixed_shift;
				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			complex operator*(const complex& other)const& noexcept{
				__int128_t Real_result=((__int128_t)Real*other.Real-(__int128_t)Imaginary*other.Imaginary)>>Fixed_shift,
						   Imaginary_result=((__int128_t)Real*other.Imaginary+(__int128_t)Imaginary*other.Real)>>Fixed_shift;
				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}


			friend complex operator/(const complex& c,const double d)noexcept{
				__int128_t Real_result=(((__int128_t)c.Real<<Fixed_shift)/FixedPoint(d)),
						   Imaginary_result=(((__int128_t)c.Imaginary<<Fixed_shift)/FixedPoint(d));
				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend complex operator/(const double& d,const complex& c)noexcept{
				__int128_t denominator=((__int128_t)c.Real*c.Real+c.Imaginary*c.Imaginary);
				__int128_t Real_result=FixedPoint(d)>c.Real?
									   ((__int128_t)FixedPoint(d)<<Fixed_shift)/denominator*c.Real:
									   ((__int128_t)c.Real<<Fixed_shift)/denominator*FixedPoint(d),

						   Imaginary_result=FixedPoint(d)>c.Imaginary?
								            -((__int128_t)FixedPoint(d)<<Fixed_shift)/denominator*c.Imaginary:
											-((__int128_t)c.Imaginary<<Fixed_shift)/denominator*FixedPoint(d);

				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend complex operator/(const complex& c,const long long d)noexcept{
				__int128_t Real_result=(((__int128_t)c.Real<<Fixed_shift)/d),
						   Imaginary_result=(((__int128_t)c.Imaginary<<Fixed_shift)/d);
				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend complex operator/(const long long& d,const complex& c)noexcept{
				__int128_t denominator=((__int128_t)c.Real*c.Real+c.Imaginary*c.Imaginary);
				__int128_t Real_result=d>c.Real?
									   ((__int128_t)d<<Fixed_shift)/denominator*c.Real:
									   ((__int128_t)c.Real<<Fixed_shift)/denominator*d,

						   Imaginary_result=d>c.Imaginary?
								            -((__int128_t)d<<Fixed_shift)/denominator*c.Imaginary:
											-((__int128_t)c.Imaginary<<Fixed_shift)/denominator*d;
				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			complex operator/(const complex& other)const& noexcept{
				__int128_t denominator=((__int128_t)other.Real*other.Real+(__int128_t)other.Imaginary*other.Imaginary);
				__int128_t Real_result=(Real>other.Real?
									   (((__int128_t)Real<<Fixed_shift)/denominator*other.Real):
									   (((__int128_t)other.Real<<Fixed_shift)/denominator*Real))
									  +(Imaginary>other.Imaginary?
									   (((__int128_t)Imaginary<<Fixed_shift)/denominator*other.Imaginary):
									   (((__int128_t)other.Imaginary<<Fixed_shift)/denominator*Imaginary)),
						   Imaginary_result=(Imaginary>other.Real?
								   	   	    (((__int128_t)Imaginary<<Fixed_shift)/denominator*other.Real):
										    (((__int128_t)other.Real<<Fixed_shift)/denominator*Imaginary))
										   -(Real>other.Imaginary?
											(((__int128_t)Real<<Fixed_shift)/denominator*other.Imaginary):
											(((__int128_t)other.Imaginary<<Fixed_shift)/denominator*Real));

				return complex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}

	};

	long long int AbsoluteValue(complex C){

		__int128_t Result=((__int128_t)C.Real*C.Real+(__int128_t)C.Imaginary*C.Imaginary);

		return NewtonSqrt(Result);

	}

	private:

	int Qubit_Amount;
	std::vector<std::pair<bool,int>> Qubit_Set_Observation;
	std::unordered_map<int,std::pair<std::vector<int>,std::vector<complex>>> Entangled_Qubit_Set_Pointer;
	std::unordered_map<int,std::pair<std::vector<int>,std::vector<complex>>*> Entangled_Qubit_Set;

	void CombineEntangledQubitSet(int situation1,int situation2){

		if(Entangled_Qubit_Set.find(situation1)==Entangled_Qubit_Set.end()||Entangled_Qubit_Set.find(situation2)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> First_Qubit_Name_Set=Entangled_Qubit_Set[situation1]->first,
						 Second_Qubit_Name_Set=Entangled_Qubit_Set[situation2]->first;

		for(unsigned long long count=0;count<Second_Qubit_Name_Set.size();count++)
			First_Qubit_Name_Set.push_back(Second_Qubit_Name_Set[count]);

		Entangled_Qubit_Set[situation1]->first=First_Qubit_Name_Set;

		std::vector<complex> First_Qubit_Set=Entangled_Qubit_Set[situation1]->second,
							 Second_Qubit_Set=Entangled_Qubit_Set[situation2]->second,
							 Result(First_Qubit_Set.size()*Second_Qubit_Set.size(),complex());

		for(unsigned long long x=0;x<First_Qubit_Set.size();x++)
			for(unsigned long long y=0;y<Second_Qubit_Set.size();y++)
				if((First_Qubit_Set[x].Real!=0||First_Qubit_Set[x].Imaginary!=0)&&(Second_Qubit_Set[y].Real!=0||Second_Qubit_Set[y].Imaginary!=0))
					Result[x*Second_Qubit_Set.size()+y]=First_Qubit_Set[x]*Second_Qubit_Set[y];

		Entangled_Qubit_Set[situation1]->second=Result;

		for(unsigned long long count=0;count<Second_Qubit_Name_Set.size();count++)
			if(Second_Qubit_Name_Set[count]!=situation1)
				Entangled_Qubit_Set_Pointer.erase(Second_Qubit_Name_Set[count]);

		for(unsigned long long count=0;count<Second_Qubit_Name_Set.size();count++)
			Entangled_Qubit_Set[Second_Qubit_Name_Set[count]]=Entangled_Qubit_Set[situation1];

	}
	constexpr inline unsigned long long BitRemova(unsigned long long num,const unsigned long long situation){

		return (num&((1ULL<<situation)-1ULL))|((num>>(situation+1ULL))<<situation);

	}
	constexpr inline unsigned long long BitRemova(unsigned long long num,unsigned long long situation1,unsigned long long situation2){

		if(situation1<situation2){
			situation1=situation1^situation2;
			situation2=situation1^situation2;
			situation1=situation1^situation2;
		}

		num=(num&((1ULL<<situation1)-1ULL))|((num>>(situation1+1ULL))<<situation1);

		return (num&((1ULL<<situation2)-1ULL))|((num>>(situation2+1ULL))<<situation2);

	}
	constexpr inline unsigned long long BitRemova(unsigned long long num,unsigned long long situation1,unsigned long long situation2,unsigned long long situation3){

		if(situation1<situation2){
			situation1=situation1^situation2;
			situation2=situation1^situation2;
			situation1=situation1^situation2;
		}
		if(situation1<situation3){
			situation1=situation1^situation3;
			situation3=situation1^situation3;
			situation1=situation1^situation3;
		}

		num=(num&((1ULL<<situation1)-1ULL))|((num>>(situation1+1ULL))<<situation1);

		if(situation2<situation3){
			situation2=situation2^situation3;
			situation3=situation2^situation3;
			situation2=situation2^situation3;
		}
		num=(num&((1ULL<<situation2)-1ULL))|((num>>(situation2+1ULL))<<situation2);

		return (num&((1ULL<<situation3)-1ULL))|((num>>(situation3+1ULL))<<situation3);

	}
	constexpr inline unsigned long long BitAdd(unsigned long long num,const unsigned long long situation,const unsigned long long add){

		return ((num>>situation)<<(situation+1ULL))|(add<<situation)|(num&((1ULL<<situation)-1ULL));

	}
	constexpr inline unsigned long long BitAdd(unsigned long long num,unsigned long long situation1,unsigned long long add1,unsigned long long situation2,unsigned long long add2){

		if(situation1>situation2){
			situation1=situation1^situation2;
			situation2=situation1^situation2;
			situation1=situation1^situation2;

			add1=add1^add2;
			add2=add1^add2;
			add1=add1^add2;
		}

		num=((num>>situation1)<<(situation1+1ULL))|(add1<<situation1)|(num&((1ULL<<situation1)-1ULL));

		return ((num>>situation2)<<(situation2+1ULL))|(add2<<situation2)|(num&((1ULL<<situation2)-1ULL));

	}
	constexpr inline unsigned long long BitAdd(unsigned long long num,unsigned long long situation1,unsigned long long add1,unsigned long long situation2,unsigned long long add2,unsigned long long situation3,unsigned long long add3){

		if(situation1>situation2){
			situation1=situation1^situation2;
			situation2=situation1^situation2;
			situation1=situation1^situation2;

			add1=add1^add2;
			add2=add1^add2;
			add1=add1^add2;
		}
		if(situation1>situation3){
			situation1=situation1^situation3;
			situation3=situation1^situation3;
			situation1=situation1^situation3;

			add1=add1^add3;
			add3=add1^add3;
			add1=add1^add3;
		}

		num=((num>>situation1)<<(situation1+1ULL))|(add1<<situation1)|(num&((1ULL<<situation1)-1ULL));

		if(situation2>situation3){
			situation2=situation2^situation3;
			situation3=situation2^situation3;
			situation2=situation2^situation3;

			add2=add2^add3;
			add3=add2^add3;
			add2=add2^add3;
		}

		num=((num>>situation2)<<(situation2+1ULL))|(add2<<situation2)|(num&((1ULL<<situation2)-1ULL));

		return ((num>>situation3)<<(situation3+1ULL))|(add3<<situation3)|(num&((1ULL<<situation3)-1ULL));

	}

	public:

	Qubit_Simulation():Qubit_Amount(0){Entangled_Qubit_Set.clear();Qubit_Set_Observation.clear();Entangled_Qubit_Set_Pointer.clear();};
	Qubit_Simulation(int num){

		Qubit_Amount=num;
		Entangled_Qubit_Set.clear();
		Qubit_Set_Observation.clear();
		Entangled_Qubit_Set_Pointer.clear();

		for(int count=0;count<num;count++)
			Entangled_Qubit_Set_Pointer[count]={{count},{complex(Fixed_Point),complex()}},
			Qubit_Set_Observation.push_back({0,-1}),
			Entangled_Qubit_Set[count]=&Entangled_Qubit_Set_Pointer[count];

	};
	void ResetQubitSet(){

		Qubit_Amount=0;
		Entangled_Qubit_Set.clear();
		Qubit_Set_Observation.clear();
		Entangled_Qubit_Set_Pointer.clear();

	}
	void OuputEntangledQubitSet(const int situation){

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return;

		if(Entangled_Qubit_Set[situation]==nullptr)
			return;

		std::vector<complex> Output_Qubit=Entangled_Qubit_Set[situation]->second;
		for(unsigned long long count=0;count<Output_Qubit.size();count++)
			std::cout<<FixedPointToDouble(Output_Qubit[count].Real)<<' '<<FixedPointToDouble(Output_Qubit[count].Imaginary)<<std::endl;

	}

	bool IsObservered(const int situation){

		return Qubit_Set_Observation[situation].first;

	}

	void CutQubitSet(const int situation,long long int Probability0,long long int Probability1){

		if(Probability0==0)
			Qubit_Set_Observation[situation]={1,1},Probability1=Fixed_Point;
		else if(Probability0==1)
			Qubit_Set_Observation[situation]={1,0},Probability0=Fixed_Point;
		else{

			std::random_device rd;
			std::mt19937_64 gen(rd());
			std::uniform_int_distribution<long long int> dist(0,Fixed_Point);
			long long int RandNum=dist(gen);

			if(RandNum<=Probability0)
				Qubit_Set_Observation[situation]={1,0};
			else
				Qubit_Set_Observation[situation]={1,1};

		}

		long long int State=Qubit_Set_Observation[situation].second;
		long long int Probability=(State==0)?Probability0:Probability1;

		std::vector<complex> Qubit_Set_Value=Entangled_Qubit_Set[situation]->second,
							 New_Qubit_Set_Value((Qubit_Set_Value.size()>>1)-(Qubit_Set_Value.size()==2),complex());

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
			for(unsigned long long count=0;count<(Qubit_Set.size());count++)
				if(Qubit_Set[count]==situation)
					Qubit_Situation=count;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		for(unsigned long long int count=0;count<(Qubit_Set_Value.size()>>1)-(Qubit_Set_Value.size()==2);count++)
			if(Qubit_Set_Value[BitAdd(count,Qubit_Situation,State)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Situation,State)].Imaginary!=0)
				New_Qubit_Set_Value[count]=Qubit_Set_Value[BitAdd(count,Qubit_Situation,State)]/NewtonSqrt((__int128_t)Probability<<Fixed_shift);

		Entangled_Qubit_Set[situation]->second=New_Qubit_Set_Value;

		Qubit_Set.erase(std::remove(Qubit_Set.begin(),Qubit_Set.end(),situation),Qubit_Set.end());
		Entangled_Qubit_Set[situation]->first=Qubit_Set;

		Entangled_Qubit_Set[situation]=nullptr;

	}

	int ObserverQubit(const int situation){

		if(IsObservered(situation))
			return Qubit_Set_Observation[situation].second;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return -1;

		std::pair<std::vector<int>,std::vector<complex>>* Qubit_Set_Pointer=Entangled_Qubit_Set[situation];
		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return -1;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		std::vector<complex> Qubit_Set_Value=Entangled_Qubit_Set[situation]->second;
		long long int Probability0=0,Probability1=0;

		for(unsigned long long int count=0;count<Qubit_Set_Value.size()>>1;count++){

			if(Qubit_Set_Value[BitAdd(count,Qubit_Situation,0ULL)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Situation,0ULL)].Imaginary!=0){
				__int128_t Caculate=AbsoluteValue(Qubit_Set_Value[BitAdd(count,Qubit_Situation,0ULL)]);
				Caculate*=Caculate;
				Caculate>>=Fixed_shift;
				Probability0+=Caculate;
			}

			if(Qubit_Set_Value[BitAdd(count,Qubit_Situation,1ULL)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Situation,1ULL)].Imaginary!=0){
				__int128_t Caculate=AbsoluteValue(Qubit_Set_Value[BitAdd(count,Qubit_Situation,1ULL)]);
				Caculate*=Caculate;
				Caculate>>=Fixed_shift;
				Probability1+=Caculate;
			}

		}

		CutQubitSet(situation,Probability0,Probability1);
		Qubit_Set=Qubit_Set_Pointer->first;
		Qubit_Set_Value=Qubit_Set_Pointer->second;

		for(long long int setting=Qubit_Set.size()-1;setting>=0;setting--){

			Qubit_Set=Qubit_Set_Pointer->first;
			Qubit_Set_Value=Qubit_Set_Pointer->second;

			Probability0=0,Probability1=0;

			for(unsigned long long int count=0;count<Qubit_Set_Value.size()>>1;count++){

				if(Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,0ULL)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,0ULL)].Imaginary!=0){
					__int128_t Caculate=AbsoluteValue(Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,0ULL)]);
					Caculate*=Caculate;
					Caculate>>=Fixed_shift;
					Probability0+=Caculate;
				}

				if(Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,1ULL)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,1ULL)].Imaginary!=0){
					__int128_t Caculate=AbsoluteValue(Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,1ULL)]);
					Caculate*=Caculate;
					Caculate>>=Fixed_shift;
					Probability1+=Caculate;
				}

			}

			if(Probability0==0||Probability1==0)
				CutQubitSet(Qubit_Set[setting],Probability0,Probability1);

		}

		return Qubit_Set_Observation[situation].second;

	}

	int GenerateQubit(){

		Qubit_Amount++;

		Entangled_Qubit_Set_Pointer[Qubit_Amount-1]={{Qubit_Amount-1},{complex(Fixed_Point),complex()}};
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Entangled_Qubit_Set.size()];
		Qubit_Set_Observation.push_back({0,-1});

		return Qubit_Amount-1;

	}

	int GeneratePositiveBellState(){

		Qubit_Amount+=2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount-2]={{Qubit_Amount-2,Qubit_Amount-1},{complex(Root_half),complex(),complex(),complex(Root_half)}};
		Entangled_Qubit_Set[Qubit_Amount-2]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];

		Qubit_Set_Observation.push_back({0,-1});
		Qubit_Set_Observation.push_back({0,-1});
		return Qubit_Amount-2;

	}

	int GenerateNegativeBellState(){

		Qubit_Amount+=2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount-2]={{Qubit_Amount-2,Qubit_Amount-1},{complex(Root_half),complex(),complex(),complex(-Root_half)}};
		Entangled_Qubit_Set[Qubit_Amount-2]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];

		Qubit_Set_Observation.push_back({0,-1});
		Qubit_Set_Observation.push_back({0,-1});
		return Qubit_Amount-2;

	}

	int GeneratePositiveOrthogonalBellState(){

		Qubit_Amount+=2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount-2]={{Qubit_Amount-2,Qubit_Amount-1},{complex(),complex(Root_half),complex(Root_half),complex()}};
		Entangled_Qubit_Set[Qubit_Amount-2]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];

		Qubit_Set_Observation.push_back({0,-1});
		Qubit_Set_Observation.push_back({0,-1});
		return Qubit_Amount-2;

	}

	int GenerateSingletBellState(){

		Qubit_Amount+=2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount-2]={{Qubit_Amount-2,Qubit_Amount-1},{complex(),complex(Root_half),complex(-Root_half),complex()}};
		Entangled_Qubit_Set[Qubit_Amount-2]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];

		Qubit_Set_Observation.push_back({0,-1});
		Qubit_Set_Observation.push_back({0,-1});
		return Qubit_Amount-2;

	}

	void Hadamard(const int situation){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++){

			 complex First=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)],
					 Second=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)];

			 Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)]=(First*(Root_half)+(Second*(Root_half)));
			 Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]=(First*(Root_half))-(Second*(Root_half));

		}

	}

	void PauliX(const int situation){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return ;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++)
			std::swap(Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)],
					  Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]);

	}

	void PauliY(const int situation){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return ;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++){

			std::swap(Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)],
					  Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]);

			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)]=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)]*complex(0LL,-Fixed_Point);
			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]*complex(0LL,Fixed_Point);

		}

	}

	void PauliZ(const int situation){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return ;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++)
			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]*(-Fixed_Point);

	}

	void PhaseS(const int situation){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return ;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++)
			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]*complex(0ULL,Fixed_Point);

	}

	void PhaseT(const int situation){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int  Qubit_Situation=-1;
		for(unsigned long long  count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return ;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++)
			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]*complex(Root_half,Root_half);


	}

	void PhaseRz(const int situation,const double Angle){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return ;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		double Cos_Angle=cos(Angle*Pi/180/2),
			   Sin_Angle=sin(Angle*Pi/180/2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++){

			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)]=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)]*complex(FixedPoint(Cos_Angle),FixedPoint(-Sin_Angle));
			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]*complex(FixedPoint(Cos_Angle),FixedPoint(Sin_Angle));

		}

	}

	void PhaseRx(const int situation,const double Angle){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return ;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		double Cos_Angle=cos(Angle*Pi/180/2),
			   Sin_Angle=sin(Angle*Pi/180/2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++){

			complex First=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)],
					Second=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)];

			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)]=First*FixedPoint(Cos_Angle)+Second*complex(0,FixedPoint(-Sin_Angle));
			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]=First*complex(0,FixedPoint(-Sin_Angle))+Second*FixedPoint(Cos_Angle);

		}

	}

	void PhaseRy(const int situation,const double Angle){

		if(IsObservered(situation))
			return ;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return ;

		Qubit_Situation=Qubit_Set.size()-Qubit_Situation-1;

		double Cos_Angle=cos(Angle*Pi/180/2),
			   Sin_Angle=sin(Angle*Pi/180/2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation]->second.size()>>1;count++){

			complex First=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)],
					Second=Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)];

			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,0ULL)]=First*FixedPoint(Cos_Angle)+Second*FixedPoint(-Sin_Angle);
			Entangled_Qubit_Set[situation]->second[BitAdd(count,Qubit_Situation,1ULL)]=First*FixedPoint(-Sin_Angle)+Second*FixedPoint(Cos_Angle);

		}

	}



	void CNOT(const int situation1,const int situation2){

		if(IsObservered(situation1)||IsObservered(situation2))
			return ;

		if(Entangled_Qubit_Set.find(situation1)==Entangled_Qubit_Set.end()||Entangled_Qubit_Set.find(situation2)==Entangled_Qubit_Set.end()||situation1==situation2)
			return ;

		if(Entangled_Qubit_Set[situation1]->first!=Entangled_Qubit_Set[situation2]->first)
			CombineEntangledQubitSet(situation1,situation2);

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation1]->first;

		int Qubit_Situation1=-1,Qubit_Situation2=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation1)
				Qubit_Situation1=count;
			else if(Qubit_Set[count]==situation2)
				Qubit_Situation2=count;

		if(Qubit_Situation1==-1||Qubit_Situation2==-1)
			return ;

		Qubit_Situation1=Qubit_Set.size()-Qubit_Situation1-1;
		Qubit_Situation2=Qubit_Set.size()-Qubit_Situation2-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation1]->second.size()>>2;count++)
			std::swap(Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,1ULL)],Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,0ULL)]);

	}

	void CZ(const int situation1,const int situation2){

		if(IsObservered(situation1)||IsObservered(situation2))
			return ;

		if(Entangled_Qubit_Set.find(situation1)==Entangled_Qubit_Set.end()||Entangled_Qubit_Set.find(situation2)==Entangled_Qubit_Set.end()||situation1==situation2)
			return ;

		if(Entangled_Qubit_Set[situation1]->first!=Entangled_Qubit_Set[situation2]->first)
			CombineEntangledQubitSet(situation1,situation2);

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation1]->first;

		int Qubit_Situation1=-1,Qubit_Situation2=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation1)
				Qubit_Situation1=count;
			else if(Qubit_Set[count]==situation2)
				Qubit_Situation2=count;

		if(Qubit_Situation1==-1||Qubit_Situation2==-1)
			return ;

		Qubit_Situation1=Qubit_Set.size()-Qubit_Situation1-1;
		Qubit_Situation2=Qubit_Set.size()-Qubit_Situation2-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation1]->second.size()>>2;count++)
			Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,1ULL)]=Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,1ULL)]*complex(-Fixed_Point);

	}

	void SWAP(const int situation1,const int situation2){

		if(IsObservered(situation1)||IsObservered(situation2))
			return ;

		if(Entangled_Qubit_Set.find(situation1)==Entangled_Qubit_Set.end()||Entangled_Qubit_Set.find(situation2)==Entangled_Qubit_Set.end()||situation1==situation2)
			return ;

		if(Entangled_Qubit_Set[situation1]->first!=Entangled_Qubit_Set[situation2]->first)
			CombineEntangledQubitSet(situation1,situation2);

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation1]->first;

		int Qubit_Situation1=-1,Qubit_Situation2=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation1)
				Qubit_Situation1=count;
			else if(Qubit_Set[count]==situation2)
				Qubit_Situation2=count;

		if(Qubit_Situation1==-1||Qubit_Situation2==-1)
			return ;

		Qubit_Situation1=Qubit_Set.size()-Qubit_Situation1-1;
		Qubit_Situation2=Qubit_Set.size()-Qubit_Situation2-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation1]->second.size()>>2;count++)
			std::swap(Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,0ULL)],Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,0ULL,Qubit_Situation2,1ULL)]);

	}

	void iSWAP(const int situation1,const int situation2){

		if(IsObservered(situation1)||IsObservered(situation2))
			return ;

		if(Entangled_Qubit_Set.find(situation1)==Entangled_Qubit_Set.end()||Entangled_Qubit_Set.find(situation2)==Entangled_Qubit_Set.end()||situation1==situation2)
			return ;

		if(Entangled_Qubit_Set[situation1]->first!=Entangled_Qubit_Set[situation2]->first)
			CombineEntangledQubitSet(situation1,situation2);

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation1]->first;

		int Qubit_Situation1=-1,Qubit_Situation2=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation1)
				Qubit_Situation1=count;
			else if(Qubit_Set[count]==situation2)
				Qubit_Situation2=count;

		if(Qubit_Situation1==-1||Qubit_Situation2==-1)
			return ;

		Qubit_Situation1=Qubit_Set.size()-Qubit_Situation1-1;
		Qubit_Situation2=Qubit_Set.size()-Qubit_Situation2-1;

		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation1]->second.size()>>2;count++){
			std::swap(Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,0ULL)],
					  Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,0ULL,Qubit_Situation2,1ULL)]);
			Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,0ULL)]=Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,0ULL)]*complex(0,Fixed_Point);
			Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,0ULL,Qubit_Situation2,1ULL)]=Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,0ULL,Qubit_Situation2,1ULL)]*complex(0,Fixed_Point);

		}
	}



	void CSWAP(const int situation1,const int situation2,const int situation3){

		if(IsObservered(situation1)||IsObservered(situation2)||IsObservered(situation3))
			return ;

		if(Entangled_Qubit_Set.find(situation1)==Entangled_Qubit_Set.end()||
		   Entangled_Qubit_Set.find(situation2)==Entangled_Qubit_Set.end()||
		   Entangled_Qubit_Set.find(situation3)==Entangled_Qubit_Set.end()||
		   situation1==situation2||
		   situation1==situation3||
		   situation1==situation3)
			return ;

		if(Entangled_Qubit_Set[situation1]->first!=Entangled_Qubit_Set[situation2]->first)
			CombineEntangledQubitSet(situation1,situation2);

		if(Entangled_Qubit_Set[situation1]->first!=Entangled_Qubit_Set[situation3]->first)
			CombineEntangledQubitSet(situation1,situation3);

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation1]->first;

		int Qubit_Situation1=-1,Qubit_Situation2=-1,Qubit_Situation3=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation1)
				Qubit_Situation1=count;
			else if(Qubit_Set[count]==situation2)
				Qubit_Situation2=count;
			else if(Qubit_Set[count]==situation3)
				Qubit_Situation3=count;

		if(Qubit_Situation1==-1||Qubit_Situation2==-1||Qubit_Situation3==-1)
			return ;

		Qubit_Situation1=Qubit_Set.size()-Qubit_Situation1-1;
		Qubit_Situation2=Qubit_Set.size()-Qubit_Situation2-1;
		Qubit_Situation3=Qubit_Set.size()-Qubit_Situation3-1;


		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation1]->second.size()>>3;count++)
			std::swap(Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,1ULL,Qubit_Situation3,0ULL)],
					  Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,0ULL,Qubit_Situation3,1LL)]);


	}

	void CCNOT(const int situation1,const int situation2,const int situation3){

		if(IsObservered(situation1)||IsObservered(situation2)||IsObservered(situation3))
			return ;

		if(Entangled_Qubit_Set.find(situation1)==Entangled_Qubit_Set.end()||
		   Entangled_Qubit_Set.find(situation2)==Entangled_Qubit_Set.end()||
		   Entangled_Qubit_Set.find(situation3)==Entangled_Qubit_Set.end()||
		   situation1==situation2||
		   situation1==situation3||
		   situation1==situation3)
			return ;

		if(Entangled_Qubit_Set[situation1]->first!=Entangled_Qubit_Set[situation2]->first)
			CombineEntangledQubitSet(situation1,situation2);

		if(Entangled_Qubit_Set[situation1]->first!=Entangled_Qubit_Set[situation3]->first)
			CombineEntangledQubitSet(situation1,situation3);

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation1]->first;

		int Qubit_Situation1=-1,Qubit_Situation2=-1,Qubit_Situation3=-1;
		for(unsigned long long count=0;count<(Qubit_Set.size());count++)
			if(Qubit_Set[count]==situation1)
				Qubit_Situation1=count;
			else if(Qubit_Set[count]==situation2)
				Qubit_Situation2=count;
			else if(Qubit_Set[count]==situation3)
				Qubit_Situation3=count;

		if(Qubit_Situation1==-1||Qubit_Situation2==-1||Qubit_Situation3==-1)
			return ;

		Qubit_Situation1=Qubit_Set.size()-Qubit_Situation1-1;
		Qubit_Situation2=Qubit_Set.size()-Qubit_Situation2-1;
		Qubit_Situation3=Qubit_Set.size()-Qubit_Situation3-1;


		for(unsigned long long count=0;count<Entangled_Qubit_Set[situation1]->second.size()>>3;count++)
			std::swap(Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,1ULL,Qubit_Situation3,0ULL)],
					  Entangled_Qubit_Set[situation1]->second[BitAdd(count,Qubit_Situation1,1ULL,Qubit_Situation2,1ULL,Qubit_Situation3,1LL)]);


	}


};

int main() {

	Qubit_Simulation a=Qubit_Simulation(3);


	a.Hadamard(0);
	a.CNOT(0,1);
	a.Hadamard(0);

	a.OuputEntangledQubitSet(0);

	std::cout<<std::endl;
	std::cout<<a.ObserverQubit(0)<<std::endl;
	std::cout<<a.ObserverQubit(1)<<std::endl;
/*
	a.PositiveGenerateBellState();

	a.Hadamard(1);
	a.CNOT(0,1);

	//a.CCNOT(0,1,2);
	//a.CSWAP(0,1,2);
	a.OuputEntangledQubitSet(0);*/


}
