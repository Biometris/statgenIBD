#ifndef POPNAME_HEADER
#define POPNAME_HEADER

#include <string>

class Pop
{
public:
	enum Type 
	{
		DH,
		Fx,
		FxDH,
		BCx,
		BCxDH,
		BCSx,
		BCSxDH,
		C3Sx,
		C3SxDH,
		C4Sx,
		C4SxDH
	};

	Pop(std::string& name);
	Type GetType() const { return type; }
	int get_x() const {return x; }

private:
    int x;
	Type type;
};

int just_test();

#endif
