#include <iostream>
#include <fstream>
#include <Winsock2.h>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <numeric>
using namespace std;

#pragma comment(lib, "ws2_32.lib")
struct CSVRow {
	int column1;
	double column2;
	double column3;
	int column4;
	double column5;
	double column6;
	double column7;
};

struct PacketData {
	double col1, col2, col3, col4, col5, col6, col7;
};

// 函数：根据给定的x值输出相应矩阵
std::vector<CSVRow> getMatrixForX(const std::vector<CSVRow>& data, int x) {
	std::vector<CSVRow> result;
	for (const auto& row : data) {
		if (row.column1 == x) {
			result.push_back(row);
		}
	}
	return result;
}

struct Row {
	int col1;
	double col2;
	double col3;
	double col4;
	double col5;
	double col7;
};

int main() {
	// 打开CSV文件
	std::ifstream csvFile("accumulated_data1207_1.csv");

	// 检查文件是否成功打开
	if (!csvFile.is_open()) {
		std::cerr << "Error opening the CSV file." << std::endl;
		return 1;
	}

	// 使用vector来存储CSV表格的行
	std::vector<CSVRow> csvData;

	// 读取CSV文件数据
	std::string line;
	while (std::getline(csvFile, line)) {
		std::istringstream iss(line);
		CSVRow row;

		// 使用逗号作为分隔符拆分每一行的数据
		char comma;
		iss >> row.column1 >> comma >> row.column2 >> comma >> row.column3 >> comma >> row.column4 >> comma >> row.column5 >> comma >> row.column6 >> comma >> row.column7;

		// 将行添加到vector中
		csvData.push_back(row);
	}

	// 关闭CSV文件
	csvFile.close();

	std::sort(csvData.begin(), csvData.end(), [](const CSVRow& a, const CSVRow& b) {
		if (a.column1 != b.column1) {
			return a.column1 < b.column1;
		}
		if (a.column2 != b.column2) {
			return a.column2 < b.column2; // Sorting in descending order for column -2
		}
		return a.column3 < b.column3;
		});

	// Get unique values from the first column
	std::vector<int> unique_values;
	for (const auto& row : csvData) {
		if (std::find(unique_values.begin(), unique_values.end(), row.column1) == unique_values.end()) {
			unique_values.push_back(row.column1);
		}
	}


	// Create a vector of vectors to store the split matrices
	std::vector<std::vector<CSVRow>> split_matrices(unique_values.size());

	// Split the matrices based on the first column value
	for (const auto& row : csvData) {
		auto it = std::find(unique_values.begin(), unique_values.end(), row.column1);
		if (it != unique_values.end()) {
			size_t index = std::distance(unique_values.begin(), it);
			split_matrices[index].push_back(row);
		}
	}

	// Create vectors to store the final concatenated matrix
	std::vector<CSVRow> data4;

	// Concatenate matrices based on the provided logic
	std::vector<std::vector<CSVRow>> split_matrices1(32);
	std::vector<int> split_matrices2(32);

	for (int i = 0; i < 32; ++i) {
		if (i % 2 == 0) {
			split_matrices1[i] = split_matrices[(i + 2) / 2 - 1];
		}
		else {
			split_matrices1[i] = split_matrices[(i - 1) / 2 + 16];
		}
		split_matrices2[i] = split_matrices1[i].size();
	}

	// Vertically concatenate matrices
	for (const auto& matrix : split_matrices1) {
		data4.insert(data4.end(), matrix.begin(), matrix.end());
	}


	// Calculate num vector
	std::vector<int> num(32);
	num[0] = split_matrices2[0];
	for (int k = 1; k < 32; ++k) {
		num[k] = split_matrices2[k] + num[k - 1];
	}

	WSADATA wsaData;
	if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
		std::cerr << "Failed to initialize Winsock." << std::endl;
		return 1;
	}

	SOCKET udpSocket = socket(AF_INET, SOCK_DGRAM, 0);
	if (udpSocket == INVALID_SOCKET) {
		std::cerr << "Failed to create UDP socket." << std::endl;
		WSACleanup();
		return 1;
	}

	sockaddr_in localAddr;
	localAddr.sin_family = AF_INET;
	localAddr.sin_port = htons(2368);
	localAddr.sin_addr.s_addr = INADDR_ANY;

	if (bind(udpSocket, reinterpret_cast<sockaddr*>(&localAddr), sizeof(localAddr)) == SOCKET_ERROR) {
		std::cerr << "Failed to bind UDP socket." << std::endl;
		closesocket(udpSocket);
		WSACleanup();
		return 1;
	}

	const int maxPackets = 167;
	const double A2 = 13.5;
	const double A1 = 9.6;
	int receivedPackets = 0;
	const double pi = acos(-1.0);
	std::vector<int> decimalValues(8192);
	std::vector<PacketData> accumulated_data;  // Change to vector of PacketData
	const int numRows = 12;
	const int numCols = 7;
	int bin_width = 1;
	std::vector<int> edges;
	for (int i = 0; i <= 92; i += bin_width) {
		edges.push_back(i - 46);
		// std::cout << edges[i] << std::endl;
	}
	auto startTime = std::chrono::steady_clock::now();

	// 在每个UDP包的循环迭代中创建一个临时的packet_data数组

	int packetCounter = 0; // 用于计数数据包的变量
	while (true) {
		char data[8192];
		int k = 1;
		std::vector<PacketData> packet_data(384);
		int bytesReceived = recvfrom(udpSocket, data, sizeof(data), 0, nullptr, nullptr);


		if (bytesReceived == SOCKET_ERROR) {
			std::cerr << "Failed to receive data." << std::endl;
		}
		else {
			decimalValues.clear();

			for (int i = 0; i < bytesReceived; ++i) {
				int decimalValue = static_cast<int>(static_cast<unsigned char>(data[i]));
				decimalValues.push_back(decimalValue);
			}

			// std::vector<int> subset(decimalValues.begin(), decimalValues.begin() + 1200);
			std::vector<int> subset;
			subset.reserve(1200);  // 提前分配空间，以避免不必要的重新分配

			// 然后给subset赋值
			subset.assign(decimalValues.begin(), decimalValues.begin() + 1200);




			for (int p = 1; p <= 12; ++p) {
				// Extract newstr
				std::vector<int> newstr(100);  // 初始化大小为100的newstr
				// std::cout << "p: " << p << std::endl;
				std::copy(subset.begin() + 100 * (p - 1), subset.begin() + 100 * p, newstr.begin());
				// std::vector<int> newstr(subset.begin() + 100 * (p - 1), subset.begin() + 100 * p);

				// 循环遍历 newstr 并输出每个元素的值
				/*std::cout << "newstr: ";
				for (int value : newstr) {
					std::cout << value << " ";
				}
				std::cout << std::endl;*/

				// cout << "hahah" << endl;
				// cout << "k=" << k << endl;
				// Loop over j
				for (int j = 1; j <= 32; ++j) {
					if (j % 2 == 1) {
						packet_data[k].col1 = -16 + (j - 1) / 2;
					}
					else {
						packet_data[k].col1 = (j - 2) / 2;
					}
					if (j % 4 == 1 || j % 4 == 2) {
						packet_data[k].col2 = round(((16 * 16 * newstr[3] + newstr[2] + 20.0 / 32 * j) / 100 + A2) * 10.0) / 10.0;
					}
					else {
						packet_data[k].col2 = round(((16 * 16 * newstr[3] + newstr[2] + 20.0 / 32 * j) / 100 + A1) * 10.0) / 10.0;
					}
					if (packet_data[k].col2 > 360.0) {
						packet_data[k].col2 -= 360.0;
					}
					// distance
					packet_data[k].col3 = (16 * 16 * newstr[5 + 3 * (j - 1)] + newstr[4 + 3 * (j - 1)]) * 0.4 / 100;
					// y
					packet_data[k].col6 = -packet_data[k].col3 * sin(packet_data[k].col2 * pi / 180.0) * cos(packet_data[k].col1 * pi / 180.0);  // y

					double packetDataCol6 = packet_data[k].col6;
					double packetDataCol3 = packet_data[k].col3;
					double packetDataCol1 = packet_data[k].col1;


					if (packetDataCol3 != 0 && packetDataCol6 >= -18 && packetDataCol6 <= 20) {
						double packetDataCol2 = packet_data[k].col2;

						for (const auto& row : split_matrices1[j - 1]) {
							// 检查第二列是否在指定的差值范围内
							if (row.column2 == packetDataCol2) {
								// 计算第三列和第一列的差值并输出
								double difference3 = std::abs(packetDataCol3 - row.column3);
								// std::cout << "在split_matrices1[" << j << "]中找到匹配: 第一列=" << row.column1 << ", 第二列=" << row.column2 << ", 差值=" << difference3 << std::endl;
								if (difference3 > 0.5) {
									// cout << "I=" << newstr[6 + 3 * (j - 1)] << endl;

									packet_data[k].col4 = newstr[6 + 3 * (j - 1)];
									packet_data[k].col5 = packetDataCol3 * cos(packetDataCol2 * pi / 180.0) * cos(packetDataCol1 * pi / 180.0);  // x
									packet_data[k].col7 = packetDataCol3 * sin(packetDataCol1 * pi / 180.0);  // z
									// std::cout << "k: " << k << std::endl;
									accumulated_data.push_back(packet_data[k]);
									// k++;
									// double packetDataCol2;
									// packetDataCol2 = newstr[7 + 3 * (j - 1)];
								}
							}
						}
					}
					else if (k != 1 && packetDataCol1 == 15 && packetDataCol3 == 0) {
						k--;
					}
					else {
						k = k;
					}

				}
			}
			// std::cout << "k=" << k << std::endl;
			// 将当前UDP包的packet_data添加到累积的数据结构中
			// accumulated_data.insert(accumulated_data.end(), packet_data.begin(), packet_data.begin()+ k);
			receivedPackets++;
			packetCounter++;

		}


		// 每167个数据包输出一次时间
		if (packetCounter == 167) {
			// cluster
			double a = accumulated_data.size();
			// 假设你有 lane 数组，用于存储符合条件的数据
			std::vector<std::vector<PacketData>> lane(4);
			int j[4] = { 0 };  // 用于追踪每个 lane 中的数据索引
			for (int i = 0; i < accumulated_data.size(); ++i) {
				double col6 = accumulated_data[i].col6;
				double col7 = accumulated_data[i].col7;

				if (col6 <= 0 && col6 >= -3.75 && col7 <= 1.25) {
					lane[0].push_back(accumulated_data[i]);
					j[0]++;
				}
				else if (col6 <= -3.75 && col6 >= -7.5 && col7 <= 1.25) {
					lane[1].push_back(accumulated_data[i]);
					j[1]++;
				}
				else if (col6 <= -7.5 && col6 >= -11.25 && col7 <= 1.25) {
					lane[2].push_back(accumulated_data[i]);
					j[2]++;
				}
				else if (col6 <= -11.25 && col6 >= -15 && col7 <= 1.25) {
					lane[3].push_back(accumulated_data[i]);
					j[3]++;
				}
			}


			// 判断点云数量是否为空
			for (int i = 0; i < lane.size(); ++i) {
				if (lane[i].empty()) {
					continue;
				}
				else {
					// 对每个 lane 按照第 5 列（假设是 col5）进行排序
					std::sort(lane[i].begin(), lane[i].end(), [](const PacketData& a, const PacketData& b) {
						return a.col5 < b.col5;
						});
				}
			}

			//// Output the contents of each lane
			//for (int i = 0; i < lane.size(); ++i) {
			//	std::cout << "Lane " << i << ": ";
			//	for (const auto& data : lane[i]) {
			//		std::cout << "(" << data.col1 << ", " << data.col2 << ", " << data.col3 << ", " << data.col4 << ", " << data.col5 << ", " << data.col6 << ", " << data.col7 << ") ";
			//	}
			//	std::cout << std::endl;
			//}


			for (int i = 0; i < lane.size(); ++i) {
				if (lane[i].empty()) {
					continue;
				}
				else {
					// 创建一个向量 counts，用于存储在每个区间内的数据点数量
					std::vector<int> counts(edges.size() - 1, 0);

					// 创建一个向量 sum_n，用于存储在每个区间前的数据点总数
					std::vector<int> sum_n;


					// 遍历每个区间，计算每个区间内的数据点数量
					for (int i_2 = 0; i_2 < edges.size() - 1; ++i_2) {
						// 如果 counts[i_2] 不为 0，说明该区间已经计算过，直接跳过
						if (counts[i_2] != 0) {
							continue;
						}

						// 使用 find_if 查找满足条件的数据点范围，返回迭代器,group_start 迭代器指向 lane[i] 中第一个满足条件
						//  packet.col5 >= edges[i_2] 的元素。同样的方式，group_end 迭代器指向 lane[i] 中第一个满足条件 packet.col5 >= edges[i_2 + 1] 的元素
						auto group_start = std::find_if(lane[i].begin(), lane[i].end(), [&edges, i_2](const PacketData& packet) {
							return packet.col5 >= edges[i_2];
							});

						auto group_end = std::find_if(lane[i].begin(), lane[i].end(), [&edges, i_2](const PacketData& packet) {
							return packet.col5 >= edges[i_2 + 1];
							});

						// 计算该区间内的数据点数量，并存储到 counts[i_2] 中
						counts[i_2] = std::distance(group_start, group_end);
						// std::cout << "counts[" << i << "] = " << i << std::endl;
						// std::cout << "counts[" << i_2 << "] = " << counts[i_2] << std::endl;

					}

					// 计算每个区间前的数据点总数
					for (int i_2 = 0; i_2 < edges.size() - 1; ++i_2) {
						sum_n.push_back(std::accumulate(counts.begin(), counts.begin() + i_2 + 1, 0));

					}

					/*std::cout << "sum_n: ";
					for (const auto& sum : sum_n) {
						std::cout << sum << " ";
					}
					std::cout << std::endl;*/

					std::vector<int> firstNonZeroIndices;
					std::vector<std::vector<int>> nonZeroGroups;
					std::vector<int> accumulatedFirstIndices;  // Store accumulated first indices

					int group_number = 0;  // 用于保存最后一个非零组的第一个值

					for (int i = 0; i < counts.size(); ++i) {
						if (counts[i] != 0) {
							firstNonZeroIndices.push_back(i);
						}
						else if (!firstNonZeroIndices.empty()) {
							// 输出每一个非零组的第一个值
							int firstIndex = firstNonZeroIndices[0];
							// Save the non-zero group's firstIndex
							accumulatedFirstIndices.push_back(firstIndex);

							// Save the non-zero group indices

							// std::cout << "非零组的第一个索引是：" << firstIndex << std::endl;
							// std::cout << "非零组的第一个值是：" << sum_n[firstIndex] << std::endl;


							group_number++;


							firstNonZeroIndices.clear();  // 清空以便下一个非零组
						}
						else if (i == counts.size() - 1) {
							accumulatedFirstIndices.push_back(91);
						}
					}

					/*std::cout << "Accumulated first indices:" << std::endl;
					for (const auto& index : accumulatedFirstIndices) {
						std::cout << index << " ";
						std::cout << "groupnumber:" << group_number << " ";
					}
					std::cout << std::endl;*/

					/*std::cout << "group_number=" << group_number << std::endl;
					std::cout << "accumulatedFirstIndices[0]=" << accumulatedFirstIndices[0] << std::endl;
					std::cout << "accumulatedFirstIndices[1]=" << accumulatedFirstIndices[1] << std::endl;*/

					// 提取 lane{i}(3:7,:) 的数据
					std::vector<PacketData> extractedData;
					for (int mm = 1; mm <= group_number - 1; mm++) {
						int startIdx = sum_n[accumulatedFirstIndices[mm - 1] - 1] + 1;
						int endIdx = sum_n[accumulatedFirstIndices[mm] - 1];

						/*std::cout << "startIdx =" << startIdx << std::endl;
						std::cout << "endIdx = " << endIdx << std::endl;*/

						// 确保 startIdx 和 endIdx 在合理范围内
						if (startIdx >= 0 && endIdx < lane[i].size()) {
							std::vector<PacketData> group(lane[i].begin() + startIdx - 1, lane[i].begin() + endIdx);

							//// 遍历 group 中的每个 PacketData，并输出其属性
							//for (const auto& packet : group) {
							//	std::cout << "col1: " << packet.col1 << ", "
							//		<< "col2: " << packet.col2 << ", "
							//		<< "col3: " << packet.col3 << ", "
							//		<< "col4: " << packet.col4 << ", "
							//		<< "col5: " << packet.col5 << ", "
							//		<< "col6: " << packet.col6 << ", "
							//		<< "col7: " << packet.col7 << std::endl;
							//}

							// Create the Excel filename!!!
							//std::string filename = "F:/SDUITC/code1/Project4/Project_cluster/object/object_" + std::to_string(i) +
							//	"_" + std::to_string(mm) + ".csv";

							//// Open the file for writing
							//std::ofstream csvFile(filename);

							//if (csvFile.is_open()) {
							//	for (const auto& packet : group) {
							//		csvFile << packet.col1 << "," << packet.col2 << "," << packet.col3 << "," << packet.col4 << "," << packet.col5 << "," << packet.col6 << "," << packet.col7 << std::endl;
							//	}
							//	csvFile.close();
							//}
						}


					}

					int lastbeginIdx = sum_n[accumulatedFirstIndices[group_number - 1] - 1] + 1;
					int lastendIdx = sum_n[accumulatedFirstIndices[group_number]];
					std::vector<PacketData> group1(lane[i].begin() + lastbeginIdx - 1, lane[i].begin() + lastendIdx);

					//// Create the Excel filename!!!
					//std::string filename = "F:/SDUITC/code1/Project4/Project_cluster/object/object_" + std::to_string(i) +
					//	"_" + std::to_string(group_number) + ".csv";

					//// Open the file for writing
					//std::ofstream csvFile(filename);

					//if (csvFile.is_open()) {
					//	for (const auto& packet : group1) {
					//		csvFile << packet.col1 << "," << packet.col2 << "," << packet.col3 << "," << packet.col4 << "," << packet.col5 << "," << packet.col6 << "," << packet.col7 << std::endl;
					//	}
					//	csvFile.close();
					//}			
				}
			}

			accumulated_data.clear();
			packet_data.clear();

			auto currentTime = std::chrono::steady_clock::now();
			auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count();

			// 计算每一帧的时间
			std::cout << elapsedTime << std::endl;

			// 重置计数器和起始时间
			packetCounter = 0;
			startTime = std::chrono::steady_clock::now();
			
		}
	}

	/*auto endTime = std::chrono::steady_clock::now();
	auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

	std::cout << "Received " << receivedPackets << " packets in " << elapsedTime << " milliseconds." << std::endl;*/

	// Open a CSV file for writing
	std::ofstream outputFile("accumulated_data1208_1.csv");

	for (const auto& data : accumulated_data) {
		outputFile << data.col1 << "," << data.col2 << "," << data.col3 << "," << data.col4 << "," << data.col5 << "," << data.col6 << "," << data.col7 << std::endl;
	}

	outputFile.close();

	std::cout.flush();

	closesocket(udpSocket);
	WSACleanup();

	return 0;
}
