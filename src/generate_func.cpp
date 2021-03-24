#include <iostream>
using namespace std;

int main(){
    cout << "#include <vector>" <<endl;
    cout << "using namespace std;"<<endl;
    cout << "\n\n"<<endl;
    for(int i = 0; i < 256; i++){
        // cout << "void func" << i << "( vector<int> &seed, int i, int j ){" <<endl;
        cout << "inline void func" << i << "( int* seed, int i, int j ){" <<endl;
        if (i){
            cout << "   int p = (i << 6) - (i << 1) + (j << 3); // p = 62*i + j*8" <<endl;
        }
        
        int cnt = 7;
        int tmp = i;
        int pos[8] = {0};
        while (tmp){
            if(tmp%2 == 1){
                pos[cnt] = 1;
            }
            cnt--;
            tmp /= 2;
        }
        for (int a = 0; a < 8; a++){
            if(pos[a] == 1)
                // cout << "   seed.push_back(p+"<<a<<");"<<endl;
                cout << "   seed[p+"<<a<<"] = 1;"<<endl;
        }
        
        // cout << "   seed.push_back(p+1);"<<endl;
        // cout << "   seed.push_back(p+2);"<<endl;
        // cout << "   seed.push_back(p+3);"<<endl;
        // cout << "   seed.push_back(p+4);"<<endl;
        // cout << "   seed.push_back(p+5);"<<endl;
        // cout << "   seed.push_back(p+6);"<<endl;
        // cout << "   seed.push_back(p+7);"<<endl;
        cout << "}"<<endl;        
    }


        cout << "\n\n"<<endl;
        // cout << "void (*find_seed_pos [])(vector<int> &seed,int i,int j) = {"<<endl;
        cout << "void (*find_seed_pos [])( int *seed,int i,int j) = {"<<endl;
        cout << "   ";
        for(int j = 0; j < 255; j++){
            cout <<"func"<<j<<",";
            if((j+1)%16 == 0){
                cout<<endl;
                cout << "   ";
            }
        }
        cout << "func255"<<endl;
        cout << "};"<<endl;
    return 0;
}