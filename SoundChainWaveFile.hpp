#pragma once

#include "SoundChain.hpp"

// TODO: Set file path
// Write more sample rates?
// Sequel with threaded write?
class WavWriterSoundChain : public SoundChainBase {
public:
	WavWriterSoundChain() {}
	~WavWriterSoundChain() {
		if (_isRecording) {
			StopRecording();
		}
	}

	void StartRecording() {
		if (_isRecording) return;

		SoundChainSettings settings = ReadSettings();

		// Open file
		audioFile.open("output.wav", std::ios::binary);

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

		_isRecording = true;
	}
	void StopRecording() {
		if (!_isRecording) return;

		int postAudioPosition = audioFile.tellp();

		audioFile.seekp(preAudioPosition - 4);
		WriteIntToFile(postAudioPosition - preAudioPosition, 4);

		audioFile.seekp(4, std::ios::beg);
		WriteIntToFile(postAudioPosition - 8, 4);

		audioFile.close();

		_isRecording = false;
	}

private:
	bool _isRecording = false;
	std::ofstream audioFile;
	int bitDepth = 16;
	int bitDepthScale = 32767;
	int preAudioPosition;

	void Reset() override {
		if (_isRecording) {
			StopRecording();
			StartRecording();
		}
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		if (_isRecording) {
			// Write sample to file
			int numberOfSamples = numberOfFrames * ReadSettings().Channels;
			for (int sample = 0; sample < numberOfSamples; sample++) {
				signed long int value = buffPtr[sample] * bitDepthScale;
				WriteIntToFile(value, bitDepth/8);
			}
		}
	}

	void WriteIntToFile(int value, int byteSize) {
		audioFile.write(reinterpret_cast<const char*>(&value), byteSize);
	}
};

// TODO: Set file path
// Read more samplerates? Currently only supports 16bit
// Sequel with threaded open?
class WavReaderSoundChain : public SoundChainBase {
public:
	WavReaderSoundChain() {
		OpenFile();
	}
	~WavReaderSoundChain() {
		if (data != nullptr) delete[] data;
	}

private:
	std::ifstream audioFile;
	float* data = nullptr;
	int dataSize = 0;
	int currentDataIndex = 0.0;
	int channels = 2;
	int sampleRate = 44100;
	int bitDepth = 16;
	float bitDepthScale = 32767.0f;

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

	void OpenFile() {
		audioFile.open("input.wav", std::ios::binary);

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

				bitDepthScale = pow(2, bitDepth) / 2 - 1;

				// std::cout << channels << " " << bitDepth << " " << bitDepthScale << std::endl;
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
		data = new float[dataSize];
		int dataIndex = 0;

		for (int i = 0; i < dataSize; i++) {
			signed short int value = 0;
			audioFile.read(reinterpret_cast<char*>(&value), bitDepth/8);

			data[dataIndex] = (float)value / bitDepthScale;
			dataIndex++;
		}
	}
};