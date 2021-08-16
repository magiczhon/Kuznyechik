#include "stdafx.h"
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<fstream>
#include <bitset>
#include<Windows.h>
#include<ctime>

using namespace std;

void FunctX(bitset<128> &vec, bitset<128> key) 
{
	/*for (int i = 0; i < 8; i++)
		vec[i] = vec[i] ^ key[i];*/
	vec = vec ^ key;
}

void Pi(bitset<8> &vec)
{
	int mas[256] = { 252, 238, 221, 17, 207, 110, 49, 22, 251, 196, 250, 218, 35, 197, 4, 77, 233, 119, 240, 219, 147, 46, 153, 186, 23, 54, 241,
		187, 20, 205, 95, 193, 249, 24, 101, 90, 226, 92, 239, 33, 129, 28, 60, 66, 139, 1, 142, 79, 5, 132, 2, 174, 227, 106, 143, 160, 6,
		11, 237, 152, 127, 212, 211, 31, 235, 52, 44, 81, 234, 200, 72, 171, 242, 42, 104, 162, 253, 58, 206, 204, 181, 112, 14, 86, 8, 12,
		118, 18, 191, 114, 19, 71, 156, 183, 93, 135, 21, 161, 150, 41, 16, 123, 154, 199, 243, 145, 120, 111, 157, 158, 178, 177, 50, 117, 25,
		61, 255, 53, 138, 126, 109, 84, 198, 128, 195, 189, 13, 87, 223, 245, 36, 169, 62, 168, 67, 201, 215, 121, 214, 246, 124, 34, 185, 3, 224,
		15, 236, 222, 122, 148, 176, 188, 220, 232, 40, 80, 78, 51, 10, 74, 167, 151, 96, 115, 30, 0, 98, 68, 26, 184, 56, 130, 100, 159, 38,
		65, 173, 69, 70, 146, 39, 94, 85, 47, 140, 163, 165, 125, 105, 213, 149, 59, 7, 88, 179, 64, 134, 172, 29, 247, 48, 55, 107, 228,
		136, 217, 231, 137, 225, 27, 131, 73, 76, 63, 248, 254, 141, 83, 170, 144, 202, 216, 133, 97, 32, 113, 103, 164, 45, 43, 9, 91, 203, 155,
		37, 208, 190, 229, 108, 82, 89, 166, 116, 210, 230, 244, 180, 192, 209, 102, 175, 194, 57, 75, 99, 182 };

	int tmp = vec.to_ulong();
	bitset<8> temp(mas[tmp]);
	vec = temp;
}

void FunctS(bitset<128> &vec)
{
	bitset<8> temp;
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			temp[j] = vec[8 * i + j ];
		}

		Pi(temp);

		for (int l = 0; l < 8; l++)
		{
			vec[8 * i + l] = temp[l];
		}
	}
}


bitset<8> MulPolinoms(bitset<8> polinom1, bitset<8> polinom2)
{
	bitset<8> balance(0b11000011);
	bitset<8> result(0x00);

	for (int i = 0; i < 8; i++)
	{
		if (polinom1[i] == 1)
		{
			result = result ^ polinom2;
		}

		if (polinom2[7] == 0)
		{
			polinom2 = polinom2 << (1);
		}
		else
		{
			polinom2 = polinom2 << (1) ^ balance;
		}
	}
	return result;
}

bitset<8> Functl(bitset<128> vec, bitset<8> ** mas)
{
	bitset<8> result;
	bitset<8> *temp = new bitset<8>[16];
	bitset<8> vect;

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			vect[j] = vec[8 * i + j];
		}
		temp[i] = vect;
	}

	result =  mas[temp[15].to_ulong()][148] ^ mas[temp[14].to_ulong()][32]  ^ mas[temp[13].to_ulong()][133] ^ mas[temp[12].to_ulong()][16] ^
			  mas[temp[11].to_ulong()][194] ^ mas[temp[10].to_ulong()][192] ^ mas[temp[9].to_ulong()][1]    ^ mas[temp[8].to_ulong()][251] ^
			  mas[temp[7].to_ulong()][1]    ^ mas[temp[6].to_ulong()][192]  ^ mas[temp[5].to_ulong()][194]  ^ mas[temp[4].to_ulong()][16]  ^
			  mas[temp[3].to_ulong()][133]  ^ mas[temp[2].to_ulong()][32]   ^ mas[temp[1].to_ulong()][148]  ^ mas[temp[0].to_ulong()][1];

	return result;
}

bitset<128> FunctR(bitset<128> vec, bitset<8> ** mas)
{
	bitset<8> temp = Functl(vec, mas);
	vec = vec >> (8);
	for (int j = 0; j < 8; j++)
		vec[120 + j] = temp[j];
	return vec;
}

void FunctL(bitset<128> &vec, bitset<8> ** mas)
{
	for (int i = 0; i < 16; i++)
		vec = FunctR(vec, mas);
}

void FunctLwithMatrixL(bitset<8> *vec, bitset<8> ** MatrixL, bitset<8> ** mas)
{
	bitset<8> * temp = new bitset<8>[16];

	for (int i = 0; i < 16; i++)
		temp[i] = vec[i];

	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			vec[15 - i] = vec[15 - i] ^ MulPolinoms(MatrixL[j][i], temp[15 - j]);
		}
}

bitset<256> FunctF(bitset<256> vec, bitset<128> key, bitset<8> ** mas)
{
	bitset<128> temp1;
	bitset<128> temp2;

	for (int i = 0; i < 128; i++)
	{
		temp1[i] = vec[128 + i];
		temp2[i] = vec[i];
	}

	FunctX(temp1, key);
	FunctS(temp1);
	FunctL(temp1, mas);
	vec = vec >> (128);

	for (int i = 0; i < 128; i++)
	{
		vec[128 + i] = temp1[i] ^ temp2[i];
	}

	return vec;
}

bitset<128> FunctEncryption(bitset<128> text, bitset<256> key, bitset<8> ** mass)
{
	bitset<128> *temp = new bitset<128>[32];

	for (int i = 0; i < 32; i++)
	{
		bitset<128> vec(i + 1);
		FunctL(vec, mass);
		temp[i] = vec;
	}

	bitset<256> *mas = new bitset<256>[5];
	mas[0] = key;

	for (int j = 0; j < 4; j++)
	{
		mas[j + 1] = mas[j];
		for (int l = 0; l < 8; l++)
		{
			mas[j + 1] = FunctF(mas[j + 1], temp[j * 8 + l], mass);
		}
	}

	bitset<128> *arr = new bitset<128>[10];

	for (int k = 0; k < 5; k++)
	{
		for (int t = 0; t < 128; t++)
		{
			(arr[k * 2 + 1])[t] = (mas[k])[t];
			(arr[k * 2])[t] = (mas[k])[t + 128];
		}
	}

	for (int m = 0; m < 9; m++)
	{
		FunctX(text, arr[m]);
		FunctS(text);
		FunctL(text, mass);
	}

	FunctX(text, arr[9]);

	return text;
}

bitset<8> ** MatrixL(bitset<8> ** MatrixR)
{
	bitset<8> ** MatrixR2 = new bitset<8>*[16];
	for (int k = 0; k < 16; k++)
		MatrixR2[k] = new bitset<8>[16];
	
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			MatrixR2[i][j] = 0;
			for (int k = 0; k < 16; k++)
				MatrixR2[i][j] = MatrixR2[i][j] ^ MulPolinoms(MatrixR[i][k], MatrixR[k][j]);
		}
	

	bitset<8> ** MatrixR4 = new bitset<8>*[16];
	for (int k = 0; k < 16; k++)
		MatrixR4[k] = new bitset<8>[16];

	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			MatrixR4[i][j] = 0;
			for (int k = 0; k < 16; k++)
				MatrixR4[i][j] = MatrixR4[i][j] ^ MulPolinoms(MatrixR2[i][k], MatrixR2[k][j]);
		}
	

	bitset<8> ** MatrixR8 = new bitset<8>*[16];
	for (int k = 0; k < 16; k++)
		MatrixR8[k] = new bitset<8>[16];

	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			MatrixR8[i][j] = 0;
			for (int k = 0; k < 16; k++)
				MatrixR8[i][j] = MatrixR8[i][j] ^ MulPolinoms(MatrixR4[i][k], MatrixR4[k][j]);
		}
	

	bitset<8> ** MatrixR16 = new bitset<8>*[16];
	for (int k = 0; k < 16; k++)
		MatrixR16[k] = new bitset<8>[16];

	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			MatrixR16[i][j] = 0;
			for (int k = 0; k < 16; k++)
				MatrixR16[i][j] = MatrixR16[i][j] ^ MulPolinoms(MatrixR8[i][k], MatrixR8[k][j]);
		}

	return MatrixR16;
}

int main()
{
	bitset<8> ** mas = new bitset<8> * [256];
	for (int i = 0; i < 256; i++)
		mas[i] = new bitset<8>[256];
	bitset<8> final;
	for (int i = 0; i < 256; i++)
		for (int j = 0; j < 256; j++)
		{
			bitset<8> brat1(i);
			bitset<8> brat2(j);
			final = MulPolinoms(brat1, brat2);
			mas[i][j] = final;
		}
	
	bitset<128> result;

	bitset<128> vec;
	bitset <64> vec1(0xffeeddccbbaa9988);
	bitset <64> vec2(0x1122334455667700); 
	for (int i = 0; i < 64; i++)
		vec[i] = vec1[i];
	for (int i = 0; i < 64; i++)
		vec[i + 64] = vec2[i];

	bitset<256> key;
	bitset <64> key1(0x0123456789abcdef);
	bitset <64> key2(0xfedcba9876543210);
	bitset <64> key3(0x0011223344556677);
	bitset <64> key4(0x8899aabbccddeeff);
	for (int i = 0; i < 64; i++)
		key[i] = key1[i];
	for (int i = 0; i < 64; i++)
		key[i + 64] = key2[i];
	for (int i = 0; i < 64; i++)
		key[i + 128] = key3[i];
	for (int i = 0; i < 64; i++)
		key[i + 192] = key4[i];

	for (int i = 0; i < 1; i++)
		result = FunctEncryption(vec, key, mas);
	cout << result << endl << endl;

	unsigned int end_time = clock();
	cout << end_time << endl;

	bitset<8> egor(0b10010100);
	bitset<8> dima(0b10010100);

	cout << MulPolinoms(egor, dima) << endl;




	bitset<8> ** MatrixR = new bitset<8>*[16];
	for (int k = 0; k < 16; k++)
		MatrixR[k] = new bitset<8>[16];

	bitset<8> R148(148);
	bitset<8> R32(32);
	bitset<8> R133(133);
	bitset<8> R16(16);
	bitset<8> R194(194);
	bitset<8> R192(192);
	bitset<8> R251(251);
	bitset<8> R0(0);
	bitset<8> R1(1);

	for (int i = 0; i < 15; i++)
	{
		for (int j = 1; j < 16; j++)
		{
			MatrixR[i][j] = R0;
		}
		MatrixR[i][i + 1] = R1;
	}

	for (int i = 1; i < 16; i++)
		MatrixR[15][i] = R0;

	MatrixR[0][0] = R148;
	MatrixR[1][0] = R32;
	MatrixR[2][0] = R133;
	MatrixR[3][0] = R16;
	MatrixR[4][0] = R194;
	MatrixR[5][0] = R192;
	MatrixR[6][0] = R1;
	MatrixR[7][0] = R251;
	MatrixR[8][0] = R1;
	MatrixR[9][0] = R192;
	MatrixR[10][0] = R194;
	MatrixR[11][0] = R16;
	MatrixR[12][0] = R133;
	MatrixR[13][0] = R32;
	MatrixR[14][0] = R148;
	MatrixR[15][0] = R1;

	bitset<8> ** MatrixLL = new bitset<8>*[16];
	for (int k = 0; k < 16; k++)
		MatrixLL[k] = new bitset<8>[16];

	MatrixLL = MatrixL(MatrixR);

	return 0;
}


