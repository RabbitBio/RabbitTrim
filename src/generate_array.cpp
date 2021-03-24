#include <iostream>
using namespace std;
#include <cstdint>

int main(){
    cout << "#include <cstdint>\n"<<endl;
    cout << "const int16_t seed_table[9][8][256][9] = { "<<endl;
    for (int i = 0; i < 9; i++){
        cout << "   { "<<endl;
        for(int j =0 ; j < 8; j++){
            cout << "       { "<<endl;
            for (int t = 0; t < 256; t++){
                int tmp = t;
                int pos[8] = {0};
                int cnt = 0;
                int num = 0; // 多少个
                while (tmp){
                    if(tmp%2){
                        pos[cnt] = 1;
                        num++;
                    }
                    cnt++;
                    tmp  /= 2;
                }
                // cout<<"int array_"<<i<<"_"<<j<<"_"<<t<<"["<<num+1<<"]"<<" = { ";
                cout << "           { ";
                cout << num;
                for(int p = 0; p < 8;p++){
                    if(pos[p]){
                        cout<<", "<<i*62+j*8+p-2; // -2
                        num--;
                    }
                }
                if(t<255){
                    cout << " },"<<endl;
                }else{
                    cout << " }"<<endl;
                }
                
                // cout << "seed_table["<<i<<"]["<<j<<"]["<<t<<"] = array_"<<i<<"_"<<j<<"_"<<t<<";"<<endl;
            }
            if(j<7)
                cout << "       },"<<endl;
            else
                cout << "       }"<<endl;
        }
        if(i < 8)
            cout << "   },"<<endl;
        else
            cout << "   }"<<endl;
    }
    cout << "};"<<endl;
    return 0;
}


