#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

using namespace std;

int total_ticks_read = 0;
int total_ticks_processed = 0;
int total_orders_executed = 0;



string companies[30] = {"AAPL", "AMZN", "ARISTA", "ATLASSIAN", "BANKNY", "BLKROCK", "BOFA", "CISCO", "CITI", "DELTA", "EXPEDIA", "FACEBOOK", "GOOGLE", "JPMC", "LYFT", "MICROSOFT", "MORGAN", "NETFLIX", "NVIDIA", "ORACLE", "PALOALTO", "PAYPAL", "QUALCOMM", "SALESFORCE", "SLB", "SOUTHWEST", "SPOTIFY", "TESLA", "UBER", "ZOOM"};

class PriceUpdate {
	public:
		string symbol;
		float ask,bid,ask_vol,bid_vol;
		long int timestamp;
		PriceUpdate(string symbol_p, float ask_p, float bid_p, float ask_vol_p, float bid_vol_p, const char * time) {
			symbol = symbol_p;
			ask = ask_p;
			bid = bid_p;
			ask_vol = ask_vol_p;
			bid_vol = bid_vol_p;
			int year,month,day,hour,minute,second,milisecond;
			sscanf(time, "%4d-%2d-%2d %2d:%2d:%2d.%6d", &year, &month, &day, &hour, &minute, &second, &milisecond);
			milisecond /= 1000;
			timestamp = (long int) milisecond + (second*1000) + (minute*60*1000) + (hour*60*60*1000) + (day*24*60*60*1000);				
		}
};

class SymbolData {
	public:
		string symbol;
		vector<float> ask,bid,ask_volume,bid_volume,gain,loss;
		int ups, downs, updates;
		float gain_sum,loss_sum, avg_gain, avg_loss;
		SymbolData() {}
		SymbolData(string symbol_p) {
			symbol = symbol_p;
			ups = 0;
			downs = 0;
			gain_sum = 0;
			loss_sum = 0;
			updates = 0;
			avg_gain = 0;
			avg_loss = 0;
		}
		SymbolData(SymbolData &t) {
			symbol = t.symbol;
			ask = t.ask;
			bid = t.bid;
			ask_volume = t.ask_volume;
			bid_volume = t.bid_volume;
			gain = t.gain;
			loss = t.loss;
			ups = t.ups;
			downs = t.downs;
			updates = t.updates;
			gain_sum = t.gain_sum;
			loss_sum = t.loss_sum;
			avg_gain = t.avg_gain;
			avg_loss = t.avg_loss;
		}
		void addUpdatePrice(PriceUpdate tick) {
			ask.push_back(tick.ask);
			bid.push_back(tick.bid);
			ask_volume.push_back(tick.ask_vol);
			bid_volume.push_back(tick.bid_vol);
			if(updates == 0) {
				gain.push_back(0);
				loss.push_back(0);	
			} else {
				if(bid[updates]-bid[updates-1]>0) {
					gain.push_back(bid[updates]-bid[updates-1]);
					loss.push_back(0);
					ups++;
					gain_sum += bid[updates]-bid[updates-1];
				} else if (bid[updates]-bid[updates-1]<0) {
					loss.push_back(bid[updates]-bid[updates-1]);
					gain.push_back(0);
					downs++;
					loss_sum += bid[updates]-bid[updates-1];
				} else if (ask[updates]-ask[updates-1]>0) {
					gain.push_back(ask[updates]-ask[updates-1]);
					loss.push_back(0);
					ups++;
					gain_sum += ask[updates]-ask[updates-1];
				} else {
					loss.push_back(ask[updates]-ask[updates-1]);
					gain.push_back(0);
					downs++;
					loss_sum += ask[updates]-ask[updates-1];
				}
			}
			if(updates>50) {
				if(gain[updates-50]!=0) {
					gain_sum -= gain[updates-50];
					ups--;
				} else {
					loss_sum -= loss[updates-50];
					downs--;
				}	
			}
			avg_gain = gain_sum/ups;
			avg_loss = loss_sum/downs;
			++updates;
		}
};

map<string, int> allStocksIndex;
vector<SymbolData*> allStocksData;

int main() {
    // Specify the name of your CSV file
    string filename = "combined-13.11.2023.csv";

    // Open the CSV file
    ifstream file(filename);

    // Check if the file is opened successfully
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }

    // Create empty stock data lists;
    for(int i=0;i<30;i++) {
	allStocksIndex.insert({companies[i],i});
	allStocksData.push_back(new SymbolData(companies[i]));
    }

    // Read the file line by line
    string line;
    getline(file, line);
    while (getline(file, line)) {
        // Use a string stream to parse each line
        istringstream iss(line);
        string value;
        vector<std::string> row;

        // Split the line into individual values
        while (getline(iss, value, ',')) {
            row.push_back(value.c_str());
        }

        // Use the row data
        PriceUpdate price = PriceUpdate(row[4], strtof(row[0].c_str(),NULL), strtof(row[1].c_str(), NULL), strtof(row[2].c_str(), NULL), strtof(row[3].c_str(), NULL), row[5].c_str());

	allStocksData[allStocksIndex[price.symbol]]->addUpdatePrice(price);
	total_ticks_read++;
    }

    // Close the file
    file.close();
    
    cout<<"Total ticks read "<<total_ticks_read<<endl;
    return 0;
}
