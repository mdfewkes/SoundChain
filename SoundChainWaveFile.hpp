#pragma once

#include "SoundChain.hpp"
#include <string>
#include <fstream>
#include <atomic>
#include <thread>

template <typename T, int S>
class RingBuffer {
public:
	static_assert(S > 0, "Size must be positive");

	bool Push(const T& item) {
		int headIndex = head.load(std::memory_order_relaxed);
		int next = NextIndex(headIndex);
		if (next == tail.load(std::memory_order_acquire)) {
			return false; // full
		}
		buffer[headIndex] = item;
		head.store(next, std::memory_order_release);
		return true;
	}

	T Pop() {
		T item;

		int tailIndex = tail.load(std::memory_order_relaxed);
		if (tailIndex == head.load(std::memory_order_acquire)) {
			return item; // empty
		}
		item = buffer[tailIndex];
		tail.store(NextIndex(tailIndex), std::memory_order_release);
		return item;
	}

	bool isEmpty() const {
		return head.load(std::memory_order_acquire) ==
			   tail.load(std::memory_order_acquire);
	}

private:
	constexpr int NextIndex(int i) const noexcept {
		return (i + 1) % S;
	}

	std::array<T, S+1> buffer;
	std::atomic<int> head{0};
	std::atomic<int> tail{0};
};

// TODO:
// Read more samplerates? Currently only supports 16bit
// Sequel with threaded write?
class WavWriterSoundChain : public SoundChainBase {
public:
	WavWriterSoundChain() {}
	~WavWriterSoundChain() {
		if (_isRecording.load()) {
			StopRecording();
		}
	}

	void StartRecording(std::string path = "output.wav") {
		if (_isRecording.load()) return;

		SoundChainSettings settings = ReadSettings();

		// Open file
		audioFile.open(path, std::ios::binary);

		// Write header
		audioFile << "RIFF";
		audioFile << "----";
		audioFile << "WAVE";

		// Format chunk
		audioFile <<"fmt ";
		WriteIntToFile(16, 4); //format chunk size
		WriteIntToFile(1, 2); // compression code
		WriteIntToFile(settings.Channels, 2); // number of channels
		WriteIntToFile(settings.SampleRate, 4); // Sample rate
		WriteIntToFile(settings.SampleRate * bitDepth / 8, 4); // Byte rate
		WriteIntToFile(bitDepth / 8, 2); // Block align
		WriteIntToFile(bitDepth, 2); // BitDepth

		// Data chunk
		audioFile << "data";
		audioFile << "----";

		preAudioPosition = audioFile.tellp();

		_isRecording.store(true);

		_writerThread = std::thread(WriteThread, std::ref(audioFile), std::ref(_ringBuffer), std::ref(_isRecording), std::ref(_dataReady), bitDepth);
	}
	void StopRecording() {
		if (!_isRecording.load()) return;

		_isRecording.store(false);
		_dataReady.fetch_add(1);
		_dataReady.notify_one();
		_writerThread.join();

		int postAudioPosition = audioFile.tellp();

		audioFile.seekp(preAudioPosition - 4);
		WriteIntToFile(postAudioPosition - preAudioPosition, 4);

		audioFile.seekp(4, std::ios::beg);
		WriteIntToFile(postAudioPosition - 8, 4);

		audioFile.close();
	}

	static void WriteThread(std::ofstream &audioFile, RingBuffer<float, 524288> &ringBuffer, std::atomic<bool> &isRecording, std::atomic<size_t> &dataReady, int bitDepth) {
		int lastSeenDataFlag = dataReady.load();
		float bitDepthScale = pow(2, bitDepth) / 2 - 1;
		while (isRecording.load()) {
			dataReady.wait(lastSeenDataFlag);

			while (!ringBuffer.isEmpty()) {
				float value = ringBuffer.Pop();
				switch (bitDepth){
					case 16:
						signed long int scaledValue;
						scaledValue = value * bitDepthScale;
						audioFile.write(reinterpret_cast<const char*>(&scaledValue), bitDepth/8);
						break;
					default:
						audioFile.write(0, bitDepth/8);
				}
			}

			lastSeenDataFlag = dataReady.load();
		}
	}

private:
	std::atomic<bool> _isRecording{false};
	std::atomic<size_t> _dataReady{0};
	std::ofstream audioFile;
	std::thread _writerThread;
	RingBuffer<float, 524288> _ringBuffer;
	int bitDepth = 16;
	float bitDepthScale = pow(2, bitDepth) / 2 - 1;
	int preAudioPosition;

	void Reset() override {
		if (_isRecording.load()) {
			StopRecording();
			StartRecording();
		}
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		if (_isRecording.load()) {
			// Write sample to file
			int numberOfSamples = numberOfFrames * ReadSettings().Channels;
			for (int sample = 0; sample < numberOfSamples; sample++) {
				// WriteSampleToFile(buffPtr[sample]);
				_ringBuffer.Push(sample);
			}

			_dataReady.fetch_add(1);
			_dataReady.notify_one();
		}
	}

	void WriteIntToFile(int value, int byteSize) {
		audioFile.write(reinterpret_cast<const char*>(&value), byteSize);
	}

	void WriteSampleToFile(float value) {
		switch (bitDepth){
			case 16:
				signed long int scaledValue;
				scaledValue = value * bitDepthScale;
				audioFile.write(reinterpret_cast<const char*>(&scaledValue), bitDepth/8);
				break;
			default:
				audioFile.write(0, bitDepth/8);
		}
	}
};

// TODO:
// Sequel with threaded open?
class WavReaderSoundChain : public SoundChainBase {
public:
	WavReaderSoundChain() {
		// OpenFile();
	}
	~WavReaderSoundChain() {
		CloseFile();
	}

	void OpenFile(std::string path = "input.wav") {
		if (data != nullptr) return;

		audioFile.open(path, std::ios::binary);

		while (!audioFile.eof() && data == nullptr) {
			std::string chunckID = ReadChunckID();

			if (chunckID == "RIFF") {
				int chunckSize;
				audioFile.read(reinterpret_cast<char*>(&chunckSize), 4);
				audioFile.seekg(4, std::ios::cur);
			} else if (chunckID == "fmt ") {
				int chunckSize;
				int dummy;
				audioFile.read(reinterpret_cast<char*>(&chunckSize), 4); // subchunck size
				audioFile.read(reinterpret_cast<char*>(&dummy), 2); // audio format
				audioFile.read(reinterpret_cast<char*>(&channels), 2); // num channels
				audioFile.read(reinterpret_cast<char*>(&sampleRate), 4); // sample rate
				audioFile.read(reinterpret_cast<char*>(&dummy), 4); // byte rate
				audioFile.read(reinterpret_cast<char*>(&dummy), 2); // block align
				audioFile.read(reinterpret_cast<char*>(&bitDepth), 2); // bits per sample
			} else if (chunckID == "data") {
				int chunckSize;
				audioFile.read(reinterpret_cast<char*>(&chunckSize), 4);
				ReadData(chunckSize);
				break;
			} else {
				int chunckSize;
				audioFile.read(reinterpret_cast<char*>(&chunckSize), 4);
				audioFile.seekg(chunckSize, std::ios::cur);
			}
		}

		audioFile.close();
	}

	void CloseFile() {
		if (data != nullptr) delete[] data;
	}

private:
	std::ifstream audioFile;
	float* data = nullptr;
	int dataSize = 0;
	int currentDataIndex = 0.0;
	int channels = 1;
	int sampleRate = 44100;
	int bitDepth = 16;

	void Reset() override {
		currentDataIndex = 0;
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		if (data == nullptr) return;

		int numberOfChannels = ReadSettings().Channels;
		int numberOfSamples = numberOfFrames * numberOfChannels;
		int playbackChannels = std::min(channels, numberOfChannels);

		for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
			for (int channel = 0; channel < playbackChannels; channel++) {
				if (currentDataIndex+channel >= dataSize) return;

				buffPtr[sample+channel] += data[currentDataIndex+channel];
			}
			currentDataIndex += channels;
		}
	}

	std::string ReadChunckID() {
		std::string chunckID = "";
		for (int i = 0; i < 4; i++) {
			char chunkIDChar;
			audioFile.read(reinterpret_cast<char*>(&chunkIDChar), 1);
			chunckID += chunkIDChar;
		}

		return chunckID;
	}

	void ReadData(int size) {
		dataSize = size / (bitDepth/8);
		if (data != nullptr) delete[] data;
		data = new float[dataSize];

		switch(bitDepth) {
		case 16: {
			signed short int value = 0;
			for (int i = 0; i < dataSize; i++) {
				// printf("%d / %d\n", i+1, dataSize);
				audioFile.read(reinterpret_cast<char*>(&value), 2);

				data[i] = (float)value / 32767.0f;
			}
			break;
		}
		case 24: {
			signed long int value = 0;
			unsigned char b[3];
			for (int i = 0; i < dataSize; i++) {
				audioFile.read((char*)b, 3);
				value = 
					(b[2] << 24) |
					(b[1] << 16) |
					(b[0] <<  8);

				data[i] = (float)value / 2147483647.0f;
			}
			break;
		}
		case 32: {
			signed long int value = 0;
			unsigned char b[4];
			for (int i = 0; i < dataSize; i++) {
				audioFile.read((char*)b, 4);
				value = 
					(b[3] << 24) |
					(b[2] << 16) |
					(b[1] <<  8) |
					(b[0] <<  0);

				data[i] = (float)value / 2147483647.0f;
			}
			break;
		}
		default: {
			for (int i = 0; i < dataSize; i++) {
				data[i] = 0.0f;
			}
		}
		}
	}
};