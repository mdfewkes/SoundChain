#pragma once

#include "miniaudio.h"

#include "SoundChainPlatform.hpp"

class MiniAudioSCP : public SoundChainPlatform {
public:
	MiniAudioSCP() {}
	~MiniAudioSCP() {}

	static void data_callback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
		float *outp = (float *)pOutput;
		MiniAudioSCP *scPlatform = (MiniAudioSCP *)pDevice->pUserData;
		scPlatform->FillBuffer(outp, frameCount);
	};

private:
	ma_device _device;

	void Setup() override {
		ma_device_config config;

		config = ma_device_config_init(ma_device_type_playback);
		config.playback.format = ma_format_f32;
		config.playback.channels = GetSoundChainSettings().Channels;
		config.sampleRate = GetSoundChainSettings().SampleRate;
		config.dataCallback = MiniAudioSCP::data_callback;
		config.pUserData = this;

		if (ma_device_init(NULL, &config, &_device) != MA_SUCCESS) {
			exit(EXIT_FAILURE);
		}
	}

	void Start() override {
		ma_device_start(&_device);
	}

	void End() override {
		ma_device_uninit(&_device);
	}
};